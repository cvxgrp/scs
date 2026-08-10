/* Routines modified from CSparse, T. Davis et al */

#include "csparse.h"

ScsMatrix *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                           scs_int triplet) {
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  if (!A) {
    return SCS_NULL;
  } /* out of memory */
  A->m = m; /* define dimensions and nzmax */
  A->n = n;
  A->p = (scs_int *)scs_calloc((triplet ? nzmax : n + 1), sizeof(scs_int));
  A->i = (scs_int *)scs_calloc(nzmax, sizeof(scs_int));
  A->x = values ? (scs_float *)scs_calloc(nzmax, sizeof(scs_float)) : SCS_NULL;
  return (!A->p || !A->i || (values && !A->x)) ? SCS(cs_spfree)(A) : A;
}

ScsMatrix *SCS(cs_done)(ScsMatrix *C, void *w, void *x, scs_int ok) {
  scs_free(w); /* free workspace */
  scs_free(x);
  return ok ? C : SCS(cs_spfree)(C); /* return result if OK, else free it */
}

/* C = compressed-column form of a triplet matrix T */
ScsMatrix *SCS(cs_compress)(const ScsMatrix *T, scs_int nz,
                            scs_int *idx_mapping) {
  scs_int m, n, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  scs_float *Cx, *Tx;
  ScsMatrix *C;
  m = T->m;
  n = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  C = SCS(cs_spalloc)(m, n, nz, Tx != SCS_NULL, 0); /* allocate result */
  w = (scs_int *)scs_calloc(n, sizeof(scs_int));    /* get workspace */
  if (!C || !w) {
    return SCS(cs_done)(C, w, SCS_NULL, 0);
  } /* out of memory */
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (k = 0; k < nz; k++) {
    w[Tj[k]]++; /* column counts */
  }
  SCS(cumsum)(Cp, w, n); /* column pointers */
  for (k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
    if (idx_mapping) {
      idx_mapping[k] = p;
    } /* old to new indices */
    if (Cx) {
      Cx[p] = Tx[k];
    }
  }
  return SCS(cs_done)(C, w, SCS_NULL, 1); /* success; free w and return C */
}

scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n) {
  scs_int i, nz = 0;
  scs_float nz2 = 0;
  if (!p || !c) {
    return (-1);
  } /* check inputs */
  for (i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i]; /* also in scs_float to avoid scs_int overflow */
    c[i] = p[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p[n] = nz;
  return nz2; /* return sum (c [0..n-1]) */
}

ScsMatrix *SCS(cs_spfree)(ScsMatrix *A) {
  if (!A) {
    return SCS_NULL;
  } /* do nothing if A already SCS_NULL */
  scs_free(A->p);
  scs_free(A->i);
  scs_free(A->x);
  scs_free(A);
  /* free the ScsMatrix struct and return SCS_NULL */
  return (ScsMatrix *)SCS_NULL;
}

/* Build the KKT matrix */
static const ScsSocMetric *g_soc_metric = SCS_NULL;
void SCS(set_soc_metric)(const ScsSocMetric *sm) {
  g_soc_metric = sm;
}
const ScsSocMetric *SCS(get_soc_metric)(void) {
  return g_soc_metric;
}

ScsMatrix *SCS(form_kkt)(const ScsMatrix *A, const ScsMatrix *P,
                         scs_float *diag_p, const scs_float *diag_r,
                         scs_int *diag_r_idxs, scs_int *soc_idxs,
                         scs_int borders, scs_int upper) {
  /*
   * Forms column compressed KKT matrix assumes column compressed A,P matrices.
   * Only upper OR lower triangular part is stuffed, depending on `upper` flag.
   *
   * Forms upper/lower triangular part of [(R_x + P)  A'; A -R_y]
   * Shapes: P : n x n, A: m x n.
   *
   * Output: `diag_p` will contain values of P diagonal upon completion,
   * and `diag_r_idxs` will contain the indices corresponding to the entries
   * in the returned matrix corresponding to the entries of R.
   *
   */
  scs_int h, i, j, count;
  ScsMatrix *Kcsc, *K;
  scs_int n = A->n;
  scs_int m = A->m;
  scs_int Anz = A->p[n];
  scs_int Knzmax;
  scs_int *idx_mapping;
  /* (research) SOC block-metric, hybrid layout. -W = -rI + 2r e0 e0' -
   * 2r ww' is the scalar diagonal plus a rank-2 term per block. Small
   * blocks (q < SOC_BORDER_MIN_Q) are emitted as dense -W upper
   * triangles inline: a tiny clique is harmless and the KKT dimension
   * does not grow (decisive with thousands of q=3 cones). Large blocks
   * carry the rank-2 term as two bordered columns appended after the
   * (n+m) core: an e0 border with diagonal -1/(2r) (negative side) and a
   * w border with diagonal +1/(2r) (positive side); eliminating the
   * borders reproduces -W exactly, the bordered matrix stays
   * quasi-definite, and the core keeps the scalar pattern/ordering (a
   * large clique would wreck the AMD ordering and with it the no-pivot
   * LDL accuracy). Upper layout only (the registering backend is upper). */
  const ScsSocMetric *soc =
      (soc_idxs && upper) ? SCS(get_soc_metric)() : SCS_NULL;
  scs_int soc_extra = 0, soc_nnz_total = 0, nborder = 0, b, q, st, jb, c;
  scs_int *soc_of_row = SCS_NULL, *soc_voff = SCS_NULL;
  if (soc && soc->n > 0) {
    soc_of_row = (scs_int *)scs_calloc(m, sizeof(scs_int));
    /* per-block offset into the registered vals / soc_idxs arrays (the
     * two layouts are consumed at different points of the emission) */
    soc_voff = (scs_int *)scs_calloc(soc->n, sizeof(scs_int));
    if (!soc_of_row || !soc_voff) {
      scs_free(soc_of_row);
      scs_free(soc_voff);
      return SCS_NULL;
    }
    for (i = 0; i < m; ++i) {
      soc_of_row[i] = -1;
    }
    for (b = 0; b < soc->n; ++b) {
      q = soc->sizes[b];
      soc_voff[b] = soc_nnz_total;
      if (q < SOC_BORDER_MIN_Q) {
        soc_extra += q * (q + 1) / 2 - q; /* dense block replaces diag */
        soc_nnz_total += q * (q + 1) / 2;
        for (i = 0; i < q; ++i) {
          soc_of_row[soc->starts[b] + i] = b;
        }
      } else if (borders) {
        soc_extra += q + 3; /* u1, q x u2, 2 border diagonals */
        soc_nnz_total += q + 2;
        nborder += 2;
      } else {
        soc_nnz_total += q + 2; /* keep offsets stable; nothing emitted */
      }
    }
  }
  if (P) {
    /* Upper bound P + I triangular component NNZs as Pnz + n */
    Knzmax = n + m + Anz + P->p[n] + soc_extra;
  } else {
    Knzmax = n + m + Anz + soc_extra;
  }
  K = SCS(cs_spalloc)(m + n + nborder, m + n + nborder, Knzmax, 1, 1);

#if VERBOSITY > 0
  scs_printf("forming kkt\n");
#endif
  /* Here we generate a triplet matrix and then compress to CSC */
  if (!K) {
    return SCS_NULL;
  }

  count = 0; /* element counter */
  if (P) {
    /* R_x + P in top left */
    for (j = 0; j < n; j++) { /* cols */
      diag_p[j] = 0.;
      /* empty column, add diagonal  */
      if (P->p[j] == P->p[j + 1]) {
        K->i[count] = j;
        K->p[count] = j;
        K->x[count] = diag_r[j];
        diag_r_idxs[j] = count; /* store the indices where diag_r occurs */
        count++;
      }
      for (h = P->p[j]; h < P->p[j + 1]; h++) {
        i = P->i[h]; /* row */
        if (i > j) { /* only upper triangular needed */
          break;
        }
        if (upper) {
          K->i[count] = i;
          K->p[count] = j;
        } else { /* lower triangular */
          /* P is passed in upper triangular, need to flip that here */
          K->i[count] = j; /* col -> row */
          K->p[count] = i; /* row -> col */
        }
        K->x[count] = P->x[h];
        if (i == j) {
          /* P has diagonal element */
          diag_p[j] = P->x[h];
          K->x[count] += diag_r[j];
          diag_r_idxs[j] = count; /* store the indices where diag_r occurs */
        }
        count++;
        /* reached the end without adding diagonal, do it now */
        if ((i < j) && (h + 1 == P->p[j + 1] || P->i[h + 1] > j)) {
          K->i[count] = j;
          K->p[count] = j;
          K->x[count] = diag_r[j];
          diag_r_idxs[j] = count; /* store the indices where diag_r occurs */
          count++;
        }
      }
    }
  } else {
    /* R_x in top left */
    for (j = 0; j < n; j++) {
      diag_p[j] = 0.;
      K->i[count] = j;
      K->p[count] = j;
      K->x[count] = diag_r[j];
      diag_r_idxs[j] = count; /* store the indices where diag_r occurs */
      count++;
    }
  }

  /* A in bottom left or A^T top right */
  for (j = 0; j < n; j++) { /* column */
    for (h = A->p[j]; h < A->p[j + 1]; h++) {
      if (upper) {
        K->p[count] = A->i[h] + n; /* column */
        K->i[count] = j;           /*row */
      } else {                     /* lower triangular */
        K->p[count] = j;           /* column */
        K->i[count] = A->i[h] + n; /* row */
      }
      K->x[count] = A->x[h];
      count++;
    }
  }

  /* -R_y at bottom right: scalar diagonal, except small SOC blocks
   * which are emitted as dense -W upper triangles (vals order: upper
   * triangle column-major per block) */
  for (j = 0; j < m; j++) {
    if (soc_of_row && soc_of_row[j] >= 0) {
      b = soc_of_row[j];
      st = soc->starts[b];
      if (j == st) {
        q = soc->sizes[b];
        c = soc_voff[b];
        for (h = 0; h < q; ++h) {    /* local column */
          for (i = 0; i <= h; ++i) { /* local row */
            K->i[count] = n + st + i;
            K->p[count] = n + st + h;
            K->x[count] = soc->vals ? soc->vals[c]
                                    : ((i == h) ? -diag_r[n + st + i] : 0.);
            if (i == h) {
              diag_r_idxs[n + st + i] = count;
            }
            soc_idxs[c++] = count;
            count++;
          }
        }
      }
      continue; /* non-start block rows handled at block start */
    }
    K->i[count] = j + n;
    K->p[count] = j + n;
    K->x[count] = -diag_r[j + n];
    diag_r_idxs[j + n] = count; /* store the indices where diag_r occurs */
    count++;
  }
  scs_free(soc_of_row);

  /* Large-block SOC border columns; registered vals per block, in
   * emission order: [d1 = -1/(2r), w_0 .. w_{q-1}, d2 = +1/(2r)]. The
   * u1 (e0) entry is the constant 1.0 and is not registered. */
  if (soc && soc->n > 0 && borders) {
    jb = n + m; /* next free border column */
    for (b = 0; b < soc->n; ++b) {
      q = soc->sizes[b];
      st = soc->starts[b];
      if (q < SOC_BORDER_MIN_Q) {
        continue; /* dense block already emitted above */
      }
      c = soc_voff[b];
      K->i[count] = n + st;
      K->p[count] = jb;
      K->x[count] = 1.0;
      count++;
      K->i[count] = jb;
      K->p[count] = jb;
      K->x[count] = soc->vals ? soc->vals[c]
                              : -1. / (2. * diag_r[n + st]);
      soc_idxs[c++] = count;
      count++;
      for (i = 0; i < q; ++i) {
        K->i[count] = n + st + i;
        K->p[count] = jb + 1;
        /* identity boost w = e0 when no registered values yet */
        K->x[count] = soc->vals ? soc->vals[c]
                                : ((i == 0) ? 1.0 : 0.0);
        soc_idxs[c++] = count;
        count++;
      }
      K->i[count] = jb + 1;
      K->p[count] = jb + 1;
      K->x[count] = soc->vals ? soc->vals[c]
                              : 1. / (2. * diag_r[n + st]);
      soc_idxs[c++] = count;
      count++;
      jb += 2;
    }
  }
  scs_free(soc_voff);

  idx_mapping = (scs_int *)scs_calloc(count, sizeof(scs_int));
  if (!idx_mapping) {
    SCS(cs_spfree)(K);
    return SCS_NULL;
  }
  Kcsc = SCS(cs_compress)(K, count, idx_mapping);
  SCS(cs_spfree)(K);
  if (Kcsc) {
    for (i = 0; i < m + n; i++) {
      diag_r_idxs[i] = idx_mapping[diag_r_idxs[i]];
    }
    for (i = 0; i < soc_nnz_total; i++) {
      soc_idxs[i] = idx_mapping[soc_idxs[i]];
    }
  }
  scs_free(idx_mapping);
  return Kcsc;
}
