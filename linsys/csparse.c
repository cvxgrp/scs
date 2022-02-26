/* Routines modified from CSparse, T. Davis et al */

#include "csparse.h"

csc *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                     scs_int triplet) {
  csc *A = (csc *)scs_calloc(1, sizeof(csc)); /* allocate the csc struct */
  if (!A) {
    return SCS_NULL;
  }         /* out of memory */
  A->m = m; /* define dimensions and nzmax */
  A->n = n;
  A->nzmax = nzmax = MAX(nzmax, 1);
  A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
  A->p = (scs_int *)scs_calloc((triplet ? nzmax : n + 1), sizeof(scs_int));
  A->i = (scs_int *)scs_calloc(nzmax, sizeof(scs_int));
  A->x = values ? (scs_float *)scs_calloc(nzmax, sizeof(scs_float)) : SCS_NULL;
  return (!A->p || !A->i || (values && !A->x)) ? SCS(cs_spfree)(A) : A;
}

csc *SCS(cs_done)(csc *C, void *w, void *x, scs_int ok) {
  scs_free(w); /* free workspace */
  scs_free(x);
  return ok ? C : SCS(cs_spfree)(C); /* return result if OK, else free it */
}

/* C = compressed-column form of a triplet matrix T */
csc *SCS(cs_compress)(const csc *T, scs_int *idx_mapping) {
  scs_int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  scs_float *Cx, *Tx;
  csc *C;
  m = T->m;
  n = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
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

csc *SCS(cs_spfree)(csc *A) {
  if (!A) {
    return SCS_NULL;
  } /* do nothing if A already SCS_NULL */
  scs_free(A->p);
  scs_free(A->i);
  scs_free(A->x);
  scs_free(A);
  return (csc *)SCS_NULL; /* free the csc struct and return SCS_NULL */
}

/* diag_p will contain values of P diagonal upon completion */
csc *SCS(form_kkt)(const ScsMatrix *A, const ScsMatrix *P, scs_float *diag_p,
                   const scs_float *diag_r, scs_int *diag_r_idxs,
                   scs_int upper) {
  /* ONLY UPPER or LOWER TRIANGULAR PART IS STUFFED
   * forms column compressed kkt matrix
   * assumes column compressed form A matrix
   *
   * forms upper/lower triangular part of [(I + P)  A'; A -I]
   * P : n x n, A: m x n.
   */
  scs_int h, i, j, count;
  csc *Kcsc, *K;
  scs_int n = A->n;
  scs_int m = A->m;
  scs_int Anz = A->p[n];
  scs_int Knzmax;
  scs_int *idx_mapping;
  if (P) {
    /* Upper bound P + I upper triangular component as Pnz + n */
    Knzmax = n + m + Anz + P->p[n];
  } else {
    Knzmax = n + m + Anz;
  }
  K = SCS(cs_spalloc)(m + n, m + n, Knzmax, 1, 1);

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
        } else {
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
      } else {
        K->p[count] = j;           /* column */
        K->i[count] = A->i[h] + n; /* row */
      }
      K->x[count] = A->x[h];
      count++;
    }
  }

  /* -R_y at bottom right */
  for (j = 0; j < m; j++) {
    K->i[count] = j + n;
    K->p[count] = j + n;
    K->x[count] = -diag_r[j + n];
    diag_r_idxs[j + n] = count; /* store the indices where diag_r occurs */
    count++;
  }

  K->nz = count;
  idx_mapping = (scs_int *)scs_calloc(K->nz, sizeof(scs_int));
  Kcsc = SCS(cs_compress)(K, idx_mapping);
  for (i = 0; i < m + n; i++) {
    diag_r_idxs[i] = idx_mapping[diag_r_idxs[i]];
  }
  SCS(cs_spfree)(K);
  scs_free(idx_mapping);
  return Kcsc;
}
