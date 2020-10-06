#include "private.h"

/* compressed sparse column matrix */
struct SPARSE_MATRIX /* matrix in compressed-column or triplet form */
{
  scs_int nzmax; /* maximum number of entries */
  scs_int m;     /* number of rows */
  scs_int n;     /* number of columns */
  scs_int *p;    /* column pointers (size n+1) or col indices (size nzmax) */
  scs_int *i;    /* row indices, size nzmax */
  scs_float *x;  /* numerical values, size nzmax */
  scs_int nz;    /* # of entries in triplet matrix, -1 for compressed-col */
};

char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsMatrix *P,
                              const ScsSettings *stgs) {
  char *tmp = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(tmp, "lin-sys:  sparse-direct\n\t  nnz(A): %li, nnz(P): %li\n",
          (long)A->p[A->n], P ? (long)P->p[P->n] : 0l);
  return tmp;
}

char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  scs_int n = p->L->n;
  sprintf(str, "lin-sys: nnz(L): %li\n", (long)(p->L->p[n] + n));
  return str;
}

/* wrapper for free */
static void *cs_free(void *p) {
  if (p) {
    scs_free(p);
  }                /* free p if it is not already SCS_NULL */
  return SCS_NULL; /* return SCS_NULL to simplify the use of cs_free */
}

static csc *cs_spfree(csc *A) {
  if (!A) {
    return SCS_NULL;
  } /* do nothing if A already SCS_NULL */
  cs_free(A->p);
  cs_free(A->i);
  cs_free(A->x);
  return (csc *)cs_free(A); /* free the csc struct and return SCS_NULL */
}

void SCS(free_lin_sys_work)(ScsLinSysWork *p) {
  if (p) {
    cs_spfree(p->L);
    cs_spfree(p->kkt);
    scs_free(p->perm);
    scs_free(p->Dinv);
    scs_free(p->bp);
    scs_free(p->scale_idxs);
    scs_free(p->Lnz);
    scs_free(p->iwork);
    scs_free(p->etree);
    scs_free(p->D);
    scs_free(p->bwork);
    scs_free(p->fwork);
    scs_free(p);
  }
}

static csc *cs_spalloc(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                       scs_int triplet) {
  csc *A = (csc *)scs_calloc(1, sizeof(csc)); /* allocate the csc struct */
  if (!A) {
    return SCS_NULL;
  }         /* out of memory */
  A->m = m; /* define dimensions and nzmax */
  A->n = n;
  A->nzmax = nzmax = MAX(nzmax, 1);
  A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
  A->p = (scs_int *)scs_malloc((triplet ? nzmax : n + 1) * sizeof(scs_int));
  A->i = (scs_int *)scs_malloc(nzmax * sizeof(scs_int));
  A->x = values ? (scs_float *)scs_malloc(nzmax * sizeof(scs_float)) : SCS_NULL;
  return (!A->p || !A->i || (values && !A->x)) ? cs_spfree(A) : A;
}

static csc *cs_done(csc *C, void *w, void *x, scs_int ok) {
  cs_free(w); /* free workspace */
  cs_free(x);
  return ok ? C : cs_spfree(C); /* return result if OK, else free it */
}

/* C = compressed-column form of a triplet matrix T */
static csc *cs_compress(const csc *T, scs_int *idx_mapping) {
  scs_int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  scs_float *Cx, *Tx;
  csc *C;
  m = T->m;
  n = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
  C = cs_spalloc(m, n, nz, Tx != SCS_NULL, 0);   /* allocate result */
  w = (scs_int *)scs_calloc(n, sizeof(scs_int)); /* get workspace */
  if (!C || !w) {
    return cs_done(C, w, SCS_NULL, 0);
  } /* out of memory */
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (k = 0; k < nz; k++) w[Tj[k]]++; /* column counts */
  SCS(cumsum)(Cp, w, n);               /* column pointers */
  for (k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
    idx_mapping[k] = p;         /* old to new indices */
    if (Cx) {
      Cx[p] = Tx[k];
    }
  }
  return cs_done(C, w, SCS_NULL, 1); /* success; free w and return C */
}

static csc *form_kkt(const ScsMatrix *A, const ScsMatrix *P,
                     const ScsSettings *s, scs_int *scale_idxs) {
  /* ONLY UPPER TRIANGULAR PART IS STUFFED
   * forms column compressed kkt matrix
   * assumes column compressed form A matrix
   *
   * forms upper triangular part of [(I + P)  A'; A -I]
   * P : n x n, A: m x n.
   */
  scs_int i, j, k, kk;
  csc *Kcsc, *K;
  scs_int n = A->n;
  scs_int m = A->m;
  scs_int Anz = A->p[n];
  scs_int Knzmax;
  scs_int *idx_mapping;
  if (P) {
    /* Upper bound P + I upper trianglar component as Pnz + n */
    Knzmax = n + m + Anz + P->p[n];
  } else {
    Knzmax = n + m + Anz;
  }
  K = cs_spalloc(m + n, m + n, Knzmax, 1, 1);

#if EXTRA_VERBOSE > 0
  scs_printf("forming kkt\n");
#endif
  /* Here we generate a triplet matrix and then compress to CSC */
  if (!K) {
    return SCS_NULL;
  }

  kk = 0; /* element counter */
  if (P) {
    /* I + P in top left */
    for (j = 0; j < P->n; j++) { /* cols */
      /* empty column, add diagonal  */
      if (P->p[j] == P->p[j + 1]) {
        K->i[kk] = j;
        K->p[kk] = j;
        K->x[kk] = s->rho_x;
        kk++;
      }
      for (k = P->p[j]; k < P->p[j + 1]; k++) {
        i = P->i[k]; /* row */
        if (i > j) { /* only upper triangular needed */
          break;
        }
        K->i[kk] = i;
        K->p[kk] = j;
        K->x[kk] = P->x[k];
        if (i == j) {
          /* P has diagonal element */
          K->x[kk] += s->rho_x;
        }
        kk++;
        /* reached the end without adding diagonal, do it now */
        if ((i < j) && (k + 1 == P->p[j + 1] || P->i[k + 1] > j)) {
          K->i[kk] = j;
          K->p[kk] = j;
          K->x[kk] = s->rho_x;
          kk++;
        }
      }
    }
  } else {
    /* rho_x * I in top left */
    for (k = 0; k < A->n; k++) {
      K->i[kk] = k;
      K->p[kk] = k;
      K->x[kk] = s->rho_x;
      kk++;
    }
  }

  /* A^T at top right */
  for (j = 0; j < n; j++) {
    for (k = A->p[j]; k < A->p[j + 1]; k++) {
      K->p[kk] = A->i[k] + n;
      K->i[kk] = j;
      K->x[kk] = A->x[k];
      kk++;
    }
  }

  /* -scale^-1 * I at bottom right */
  for (k = 0; k < m; k++) {
    K->i[kk] = k + n;
    K->p[kk] = k + n;
    K->x[kk] = -1. / s->scale;
    scale_idxs[k] = kk; /* store the indices where scale occurs */
    kk++;
  }
  K->nz = kk;
  idx_mapping = (scs_int *)scs_malloc(K->nz * sizeof(scs_int));
  Kcsc = cs_compress(K, idx_mapping);
  for (i = 0; i < A->m; i++) {
    scale_idxs[i] = idx_mapping[scale_idxs[i]];
  }
  cs_spfree(K);
  scs_free(idx_mapping);
  return Kcsc;
}

static scs_int _ldl_init(csc *A, scs_int *P, scs_float **info) {
  *info = (scs_float *)scs_malloc(AMD_INFO * sizeof(scs_float));
  return amd_order(A->n, A->p, A->i, P, (scs_float *)SCS_NULL, *info);
}

/* call only once */
static scs_int ldl_prepare(ScsLinSysWork *p) {
  csc *kkt = p->kkt, *L = p->L;
  scs_int n = kkt->n;
  p->etree = (scs_int *)scs_malloc(n * sizeof(scs_int));
  p->Lnz = (scs_int *)scs_malloc(n * sizeof(scs_int));
  p->iwork = (scs_int *)scs_malloc(3 * n * sizeof(scs_int));
  L->p = (scs_int *)scs_malloc((1 + n) * sizeof(scs_int));
  L->nzmax = QDLDL_etree(n, kkt->p, kkt->i, p->iwork, p->Lnz, p->etree);
  if (L->nzmax < 0) {
    scs_printf("Error in elimination tree calculation.\n");
    return L->nzmax;
  }
  L->x = (scs_float *)scs_malloc(L->nzmax * sizeof(scs_float));
  L->i = (scs_int *)scs_malloc(L->nzmax * sizeof(scs_int));
  p->Dinv = (scs_float *)scs_malloc(n * sizeof(scs_float));
  p->D = (scs_float *)scs_malloc(n * sizeof(scs_float));
  p->bwork = (scs_int *)scs_malloc(n * sizeof(scs_int));
  p->fwork = (scs_float *)scs_malloc(n * sizeof(scs_float));
  return L->nzmax;
}

/* can call many times */
static scs_int ldl_factor(ScsLinSysWork *p) {
  scs_int factor_status;
  csc * kkt = p->kkt, * L = p->L;
#if EXTRA_VERBOSE > 0
  scs_printf("numeric factorization\n");
#endif
  factor_status = QDLDL_factor(kkt->n, kkt->p, kkt->i, kkt->x, L->p, L->i, L->x,
                               p->D, p->Dinv, p->Lnz, p->etree, p->bwork,
                               p->iwork, p->fwork);
#if EXTRA_VERBOSE > 0
  scs_printf("finished numeric factorization\n");
#endif
  return factor_status;
}

static void _ldl_perm(scs_int n, scs_float *x, scs_float *b, scs_int *P) {
  scs_int j;
  for (j = 0; j < n; j++) x[j] = b[P[j]];
}

static void _ldl_permt(scs_int n, scs_float *x, scs_float *b, scs_int *P) {
  scs_int j;
  for (j = 0; j < n; j++) x[P[j]] = b[j];
}

static void _ldl_solve(scs_float *b, csc *L, scs_float *Dinv, scs_int *P,
                       scs_float *bp) {
  /* solves PLDL'P' x = b for x */
  scs_int n = L->n;
  _ldl_perm(n, bp, b, P);
  QDLDL_solve(n, L->p, L->i, L->x, Dinv, bp);
  _ldl_permt(n, b, bp, P);
}

void SCS(accum_by_atrans)(const ScsMatrix *A, ScsLinSysWork *p,
                          const scs_float *x, scs_float *y) {
  SCS(_accum_by_atrans)(A->n, A->x, A->i, A->p, x, y);
}

void SCS(accum_by_a)(const ScsMatrix *A, ScsLinSysWork *p, const scs_float *x,
                     scs_float *y) {
  SCS(_accum_by_a)(A->n, A->x, A->i, A->p, x, y, 0);
}

static scs_int *cs_pinv(scs_int const *p, scs_int n) {
  scs_int k, *pinv;
  if (!p) {
    return SCS_NULL;
  } /* p = SCS_NULL denotes identity */
  pinv = (scs_int *)scs_malloc(n * sizeof(scs_int)); /* allocate result */
  if (!pinv) {
    return SCS_NULL;
  }                                       /* out of memory */
  for (k = 0; k < n; k++) pinv[p[k]] = k; /* invert the permutation */
  return pinv;                            /* return result */
}

static csc *cs_symperm(const csc *A, const scs_int *pinv, scs_int *idx_mapping,
                       scs_int values) {
  scs_int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
  scs_float *Cx, *Ax;
  csc *C;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = cs_spalloc(n, n, Ap[n], values && (Ax != SCS_NULL), 0); /* alloc result*/
  w = (scs_int *)scs_calloc(n, sizeof(scs_int)); /* get workspace */
  if (!C || !w) {
    return cs_done(C, w, SCS_NULL, 0);
  } /* out of memory */
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (j = 0; j < n; j++) /* count entries in each column of C */
  {
    j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];
      if (i > j) {
        continue;
      }                        /* skip lower triangular part of A */
      i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
      w[MAX(i2, j2)]++;        /* column count of C */
    }
  }
  SCS(cumsum)(Cp, w, n); /* compute column pointers of C */
  for (j = 0; j < n; j++) {
    j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];
      if (i > j) {
        continue;
      }                        /* skip lower triangular part of A*/
      i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
      Ci[q = w[MAX(i2, j2)]++] = MIN(i2, j2);
      if (Cx) {
        Cx[q] = Ax[p];
      }
      idx_mapping[p] = q; /* old to new indices */
    }
  }
  return cs_done(C, w, SCS_NULL, 1); /* success; free workspace, return C */
}

static csc * permute_kkt(const ScsMatrix *A, const ScsMatrix *P,
                         const ScsSettings *stgs, ScsLinSysWork *p) {
  scs_float *info;
  scs_int *Pinv, amd_status, *idx_mapping, i;
  csc *kkt_perm, *kkt = form_kkt(A, P, stgs, p->scale_idxs);
  if (!kkt) {
    return SCS_NULL;
  }
  amd_status = _ldl_init(kkt, p->perm, &info);
  if (amd_status < 0) {
    return SCS_NULL;
  }
#if EXTRA_VERBOSE > 0
  if (stgs->verbose) {
    scs_printf("Matrix factorization info:\n");
    amd_info(info);
  }
#endif
  Pinv = cs_pinv(p->perm, A->n + A->m);
  idx_mapping = (scs_int *)scs_malloc(kkt->nzmax * sizeof(scs_int));
  kkt_perm = cs_symperm(kkt, Pinv, idx_mapping, 1);
  for (i = 0; i < A->m; i++){
    p->scale_idxs[i] = idx_mapping[p->scale_idxs[i]];
  }
  cs_spfree(kkt);
  scs_free(Pinv);
  scs_free(info);
  scs_free(idx_mapping);
  return kkt_perm;
}

scs_int SCS(should_update_scale)(scs_float factor, scs_int iter) {
  /* XXX update */
  return (factor > 5 || factor < 0.2);
}

void SCS(update_linsys_scale)(const ScsMatrix *A, const ScsMatrix *P,
                              const ScsSettings *stgs, ScsLinSysWork *p) {
  scs_int i, ldl_status;
  for (i = 0; i < A->m; ++i) {
    p->kkt->x[p->scale_idxs[i]] = -1.0 / stgs->scale;
  }
  ldl_status = ldl_factor(p);
  if (ldl_status < 0) {
    scs_printf("Error in factorize\n");
    /* XXX this is broken somehow */
    // SCS(free_lin_sys_work)(p);
    return;
  }
}

ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A, const ScsMatrix *P,
                                      const ScsSettings *stgs) {
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  scs_int n_plus_m = A->n + A->m, ldl_status, ldl_prepare_status;
  p->perm = (scs_int *)scs_malloc(sizeof(scs_int) * n_plus_m);
  p->L = (csc *)scs_malloc(sizeof(csc));
  p->bp = (scs_float *)scs_malloc(n_plus_m * sizeof(scs_float));
  p->scale_idxs = (scs_int *)scs_malloc(A->m * sizeof(scs_int));
  p->L->m = n_plus_m;
  p->L->n = n_plus_m;
  p->L->nz = -1;
  p->kkt = permute_kkt(A, P, stgs, p);
  ldl_prepare_status = ldl_prepare(p);
  ldl_status = ldl_factor(p);
  if (ldl_prepare_status < 0 || ldl_status < 0) {
    scs_printf("Error in factorize\n");
    // SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  return p;
}

scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsMatrix *P,
                           const ScsSettings *stgs, ScsLinSysWork *p,
                           scs_float *b, const scs_float *s, scs_int iter) {
  /* returns solution to linear system */
  /* Ax = b with solution stored in b */
  _ldl_solve(b, p->L, p->Dinv, p->perm, p->bp);
  return 0;
}

void SCS(normalize)(ScsMatrix *A, ScsMatrix *P, const ScsCone *k,
                    ScsScaling *scal, ScsConeWork * c) {
  SCS(_normalize)(A, P, k, scal, c);
}

void SCS(un_normalize)(ScsMatrix *A, ScsMatrix *P, const ScsScaling *scal) {
  SCS(_un_normalize)(A, P, scal);
}
