#include "private.h"

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

char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsSettings *stgs) {
  char *tmp = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(tmp, "sparse-direct, nnz in A = %li", (long)A->p[A->n]);
  return tmp;
}

char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  scs_int n = p->L->n;
  sprintf(str, "\tLin-sys: nnz in L factor: %li, avg solve time: %1.2es\n",
          (long)(p->L->p[n] + n), p->total_solve_time / (info->iter + 1) / 1e3);
  p->total_solve_time = 0;
  return str;
}

/* wrapper for free */
static void *cs_free(void *p) {
  if (p) {
    scs_free(p);
  }                /* free p if it is not already SCS_NULL */
  return SCS_NULL; /* return SCS_NULL to simplify the use of cs_free */
}

static _cs *cs_spfree(_cs *A) {
  if (!A) {
    return SCS_NULL;
  } /* do nothing if A already SCS_NULL */
  cs_free(A->p);
  cs_free(A->i);
  cs_free(A->x);
  return (_cs *)cs_free(A); /* free the _cs struct and return SCS_NULL */
}

void SCS(free_lin_sys_work)(ScsLinSysWork *p) {
  if (p) {
    cs_spfree(p->L);
    scs_free(p->P);
    scs_free(p->Dinv);
    scs_free(p->bp);
    scs_free(p);
  }
}

static _cs *cs_spalloc(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                       scs_int triplet) {
  _cs *A = (_cs *)scs_calloc(1, sizeof(_cs)); /* allocate the _cs struct */
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

static _cs *cs_done(_cs *C, void *w, void *x, scs_int ok) {
  cs_free(w); /* free workspace */
  cs_free(x);
  return ok ? C : cs_spfree(C); /* return result if OK, else free it */
}

/* C = compressed-column form of a triplet matrix T */
static _cs *cs_compress(const _cs *T) {
  scs_int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  scs_float *Cx, *Tx;
  _cs *C;
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
    if (Cx) {
      Cx[p] = Tx[k];
    }
  }
  return cs_done(C, w, SCS_NULL, 1); /* success; free w and return C */
}

static _cs *form_kkt(const ScsMatrix *A, const ScsSettings *s) {
  /* ONLY UPPER TRIANGULAR PART IS STUFFED
   * forms column compressed KKT matrix
   * assumes column compressed form A matrix
   *
   * forms upper triangular part of [I A'; A -I]
   */
  scs_int j, k, kk;
  _cs *K_cs;
  /* I at top left */
  const scs_int Anz = A->p[A->n];
  const scs_int Knzmax = A->n + A->m + Anz;
  _cs *K = cs_spalloc(A->m + A->n, A->m + A->n, Knzmax, 1, 1);

#if EXTRA_VERBOSE > 0
  scs_printf("forming KKT\n");
#endif

  if (!K) {
    return SCS_NULL;
  }
  kk = 0;
  for (k = 0; k < A->n; k++) {
    K->i[kk] = k;
    K->p[kk] = k;
    K->x[kk] = s->rho_x;
    kk++;
  }
  /* A^T at top right : CCS: */
  for (j = 0; j < A->n; j++) {
    for (k = A->p[j]; k < A->p[j + 1]; k++) {
      K->p[kk] = A->i[k] + A->n;
      K->i[kk] = j;
      K->x[kk] = A->x[k];
      kk++;
    }
  }
  /* -I at bottom right */
  for (k = 0; k < A->m; k++) {
    K->i[kk] = k + A->n;
    K->p[kk] = k + A->n;
    K->x[kk] = -1;
    kk++;
  }
  /* assert kk == Knzmax */
  K->nz = Knzmax;
  K_cs = cs_compress(K);
  cs_spfree(K);
  return K_cs;
}

static scs_int _ldl_init(_cs *A, scs_int *P, scs_float **info) {
  *info = (scs_float *)scs_malloc(AMD_INFO * sizeof(scs_float));
  return amd_order(A->n, A->p, A->i, P, (scs_float *)SCS_NULL, *info);
}

static scs_int _ldl_factor(_cs *A, _cs **L, scs_float **Dinv) {
  scs_int factor_status, n = A->n;
  scs_int *etree = (scs_int *)scs_malloc(n * sizeof(scs_int));
  scs_int *Lnz = (scs_int *)scs_malloc(n * sizeof(scs_int));
  scs_int *iwork = (scs_int *)scs_malloc(3 * n * sizeof(scs_int));
  scs_float *D, *fwork;
  scs_int *bwork;
  (*L)->p = (scs_int *)scs_malloc((1 + n) * sizeof(scs_int));
  (*L)->nzmax = QDLDL_etree(n, A->p, A->i, iwork, Lnz, etree);
  if ((*L)->nzmax < 0) {
    scs_printf("Error in elimination tree calculation.\n");
    scs_free(Lnz);
    scs_free(iwork);
    scs_free(etree);
    scs_free((*L)->p);
    return (*L)->nzmax;
  }

  (*L)->x = (scs_float *)scs_malloc((*L)->nzmax * sizeof(scs_float));
  (*L)->i = (scs_int *)scs_malloc((*L)->nzmax * sizeof(scs_int));
  *Dinv = (scs_float *)scs_malloc(n * sizeof(scs_float));
  D = (scs_float *)scs_malloc(n * sizeof(scs_float));
  bwork = (scs_int *)scs_malloc(n * sizeof(scs_int));
  fwork = (scs_float *)scs_malloc(n * sizeof(scs_float));

#if EXTRA_VERBOSE > 0
  scs_printf("numeric factorization\n");
#endif
  factor_status = QDLDL_factor(n, A->p, A->i, A->x, (*L)->p, (*L)->i, (*L)->x,
                               D, *Dinv, Lnz, etree, bwork, iwork, fwork);
#if EXTRA_VERBOSE > 0
  scs_printf("finished numeric factorization\n");
#endif

  scs_free(Lnz);
  scs_free(iwork);
  scs_free(etree);
  scs_free(D);
  scs_free(bwork);
  scs_free(fwork);
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

static void _ldl_solve(scs_float *b, _cs *L, scs_float *Dinv, scs_int *P,
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
  SCS(_accum_by_a)(A->n, A->x, A->i, A->p, x, y);
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

static _cs *cs_symperm(const _cs *A, const scs_int *pinv, scs_int values) {
  scs_int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
  scs_float *Cx, *Ax;
  _cs *C;
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
    }
  }
  return cs_done(C, w, SCS_NULL, 1); /* success; free workspace, return C */
}

static scs_int factorize(const ScsMatrix *A, const ScsSettings *stgs,
                         ScsLinSysWork *p) {
  scs_float *info;
  scs_int *Pinv, amd_status, ldl_status;
  _cs *C, *K = form_kkt(A, stgs);
  if (!K) {
    return -1;
  }
  amd_status = _ldl_init(K, p->P, &info);
  if (amd_status < 0) {
    return amd_status;
  }
#if EXTRA_VERBOSE > 0
  if (stgs->verbose) {
    scs_printf("Matrix factorization info:\n");
    amd_info(info);
  }
#endif
  Pinv = cs_pinv(p->P, A->n + A->m);
  C = cs_symperm(K, Pinv, 1);
  ldl_status = _ldl_factor(C, &p->L, &p->Dinv);
  cs_spfree(C);
  cs_spfree(K);
  scs_free(Pinv);
  scs_free(info);
  return ldl_status;
}

ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A,
                                      const ScsSettings *stgs) {
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  scs_int n_plus_m = A->n + A->m;
  p->P = (scs_int *)scs_malloc(sizeof(scs_int) * n_plus_m);
  p->L = (_cs *)scs_malloc(sizeof(_cs));
  p->bp = (scs_float *)scs_malloc(n_plus_m * sizeof(scs_float));
  p->L->m = n_plus_m;
  p->L->n = n_plus_m;
  p->L->nz = -1;

  if (factorize(A, stgs, p) < 0) {
    SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  p->total_solve_time = 0.0;
  return p;
}

scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsSettings *stgs,
                           ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_int iter) {
  /* returns solution to linear system */
  /* Ax = b with solution stored in b */
  SCS(timer) linsys_timer;
  SCS(tic)(&linsys_timer);
  _ldl_solve(b, p->L, p->Dinv, p->P, p->bp);
  p->total_solve_time += SCS(tocq)(&linsys_timer);
#if EXTRA_VERBOSE > 0
  scs_printf("linsys solve time: %1.2es\n", SCS(tocq)(&linsys_timer) / 1e3);
#endif
  return 0;
}

void SCS(normalize_a)(ScsMatrix *A, const ScsSettings *stgs, const ScsCone *k,
                      ScsScaling *scal) {
  SCS(_normalize_a)(A, stgs, k, scal);
}

void SCS(un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                         const ScsScaling *scal) {
  SCS(_un_normalize_a)(A, stgs, scal);
}
