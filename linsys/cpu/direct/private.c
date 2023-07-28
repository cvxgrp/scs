#include "private.h"

const char *scs_get_lin_sys_method(void) {
  return "sparse-direct-amd-qdldl";
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    SCS(cs_spfree)(p->L);
    SCS(cs_spfree)(p->kkt);
    scs_free(p->diag_p);
    scs_free(p->perm);
    scs_free(p->Dinv);
    scs_free(p->bp);
    scs_free(p->diag_r_idxs);
    scs_free(p->Lnz);
    scs_free(p->iwork);
    scs_free(p->etree);
    scs_free(p->D);
    scs_free(p->bwork);
    scs_free(p->fwork);
    scs_free(p);
  }
}

static scs_int _ldl_init(ScsMatrix *A, scs_int *P, scs_float **info) {
  *info = (scs_float *)scs_calloc(AMD_INFO, sizeof(scs_float));
  return amd_order(A->n, A->p, A->i, P, (scs_float *)SCS_NULL, *info);
}

/* call only once */
static scs_int ldl_prepare(ScsLinSysWork *p) {
  ScsMatrix *kkt = p->kkt, *L = p->L;
  scs_int nzmax, n = kkt->n;
  p->etree = (scs_int *)scs_calloc(n, sizeof(scs_int));
  p->Lnz = (scs_int *)scs_calloc(n, sizeof(scs_int));
  p->iwork = (scs_int *)scs_calloc(3 * n, sizeof(scs_int));
  L->p = (scs_int *)scs_calloc((1 + n), sizeof(scs_int));
  nzmax = QDLDL_etree(n, kkt->p, kkt->i, p->iwork, p->Lnz, p->etree);
  if (nzmax < 0) {
    scs_printf("Error in elimination tree calculation.\n");
    if (nzmax == -1) {
      scs_printf("Matrix is not perfectly upper triangular.\n");
    } else if (nzmax == -2) {
      scs_printf("Integer overflow in L nonzero count.\n");
    }
    return nzmax;
  }

  L->x = (scs_float *)scs_calloc(nzmax, sizeof(scs_float));
  L->i = (scs_int *)scs_calloc(nzmax, sizeof(scs_int));
  p->Dinv = (scs_float *)scs_calloc(n, sizeof(scs_float));
  p->D = (scs_float *)scs_calloc(n, sizeof(scs_float));
  p->bwork = (scs_int *)scs_calloc(n, sizeof(scs_int));
  p->fwork = (scs_float *)scs_calloc(n, sizeof(scs_float));
  return nzmax;
}

/* can call many times */
static scs_int ldl_factor(ScsLinSysWork *p, scs_int num_vars) {
  scs_int factor_status;
  ScsMatrix *kkt = p->kkt, *L = p->L;
#if VERBOSITY > 0
  scs_printf("numeric factorization\n");
#endif
  factor_status =
      QDLDL_factor(kkt->n, kkt->p, kkt->i, kkt->x, L->p, L->i, L->x, p->D,
                   p->Dinv, p->Lnz, p->etree, p->bwork, p->iwork, p->fwork);
#if VERBOSITY > 0
  scs_printf("finished numeric factorization.\n");
#endif
  if (factor_status < 0) {
    scs_printf("Error in LDL factorization when computing the nonzero "
               "elements. There are zeros in the diagonal matrix.\n");
  } else if (factor_status < num_vars) {
    scs_printf("Error in LDL factorization when computing the nonzero "
               "elements. The problem seems to be non-convex.\n");
    scs_printf("factor_status: %li, num_vars: %li\n", (long)factor_status,
               (long)num_vars);
    return -1;
  }
  p->factorizations++;
  return factor_status;
}

static void _ldl_perm(scs_int n, scs_float *x, scs_float *b, scs_int *P) {
  scs_int j;
  for (j = 0; j < n; j++)
    x[j] = b[P[j]];
}

static void _ldl_permt(scs_int n, scs_float *x, scs_float *b, scs_int *P) {
  scs_int j;
  for (j = 0; j < n; j++)
    x[P[j]] = b[j];
}

static void _ldl_solve(scs_float *b, ScsMatrix *L, scs_float *Dinv, scs_int *P,
                       scs_float *bp) {
  /* solves PLDL'P' x = b for x */
  scs_int n = L->n;
  _ldl_perm(n, bp, b, P);
  QDLDL_solve(n, L->p, L->i, L->x, Dinv, bp);
  _ldl_permt(n, b, bp, P);
}

static scs_int *cs_pinv(scs_int const *p, scs_int n) {
  scs_int k, *pinv;
  if (!p) {
    return SCS_NULL;
  } /* p = SCS_NULL denotes identity */
  pinv = (scs_int *)scs_calloc(n, sizeof(scs_int)); /* allocate result */
  if (!pinv) {
    return SCS_NULL;
  } /* out of memory */
  for (k = 0; k < n; k++)
    pinv[p[k]] = k; /* invert the permutation */
  return pinv;      /* return result */
}

static ScsMatrix *cs_symperm(const ScsMatrix *A, const scs_int *pinv,
                             scs_int *idx_mapping, scs_int values) {
  scs_int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
  scs_float *Cx, *Ax;
  ScsMatrix *C;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = SCS(cs_spalloc)(n, n, Ap[n], values && (Ax != SCS_NULL),
                      0);                        /* alloc result*/
  w = (scs_int *)scs_calloc(n, sizeof(scs_int)); /* get workspace */
  if (!C || !w) {
    return SCS(cs_done)(C, w, SCS_NULL, 0);
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
  return SCS(cs_done)(C, w, SCS_NULL,
                      1); /* success; free workspace, return C */
}

static ScsMatrix *permute_kkt(const ScsMatrix *A, const ScsMatrix *P,
                              ScsLinSysWork *p, const scs_float *diag_r) {
  scs_float *info;
  scs_int *Pinv, amd_status, *idx_mapping, i, kkt_nnz;
  ScsMatrix *kkt_perm;
  ScsMatrix *kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!kkt) {
    return SCS_NULL;
  }
  kkt_nnz = kkt->p[kkt->n];
  amd_status = _ldl_init(kkt, p->perm, &info);
  if (amd_status < 0) {
    scs_printf("AMD permutatation error.\n");
    return SCS_NULL;
  }
#if VERBOSITY > 0
  scs_printf("Matrix factorization info:\n");
  amd_info(info);
#endif
  Pinv = cs_pinv(p->perm, A->n + A->m);
  idx_mapping = (scs_int *)scs_calloc(kkt_nnz, sizeof(scs_int));
  kkt_perm = cs_symperm(kkt, Pinv, idx_mapping, 1);
  for (i = 0; i < A->n + A->m; i++) {
    p->diag_r_idxs[i] = idx_mapping[p->diag_r_idxs[i]];
  }
  SCS(cs_spfree)(kkt);
  scs_free(Pinv);
  scs_free(info);
  scs_free(idx_mapping);
  return kkt_perm;
}

void scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i, ldl_status;
  for (i = 0; i < p->n; ++i) {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i) {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }
  ldl_status = ldl_factor(p, p->n);
  if (ldl_status < 0) {
    scs_printf("Error in LDL factorization when updating.\n");
    /* TODO: this is broken somehow */
    /* SCS(free_lin_sys_work)(p); */
    return;
  }
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  scs_int n_plus_m = A->n + A->m, ldl_status, ldl_prepare_status;
  p->m = A->m;
  p->n = A->n;
  p->diag_p = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->perm = (scs_int *)scs_calloc(sizeof(scs_int), n_plus_m);
  p->L = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  p->bp = (scs_float *)scs_calloc(n_plus_m, sizeof(scs_float));
  p->diag_r_idxs = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->factorizations = 0;
  p->L->m = n_plus_m;
  p->L->n = n_plus_m;
  p->kkt = permute_kkt(A, P, p, diag_r);
  ldl_prepare_status = ldl_prepare(p);
  ldl_status = ldl_factor(p, A->n);
  if (ldl_prepare_status < 0 || ldl_status < 0) {
    scs_printf("Error in LDL initial factorization.\n");
    /* TODO: this is broken somehow */
    /* SCS(free_lin_sys_work)(p); */
    return SCS_NULL;
  }
  return p;
}

scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  /* returns solution to linear system */
  /* Ax = b with solution stored in b */
  _ldl_solve(b, p->L, p->Dinv, p->perm, p->bp);
  return 0;
}
