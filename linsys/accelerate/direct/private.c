#include "private.h"

#ifdef DLONG
#error "Apple Accelerate backend does not support DLONG (64-bit integers)"
#endif

const char *scs_get_lin_sys_method(void) {
  return "sparse-direct-apple-accelerate";
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  scs_int i, n_plus_m;
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p) {
    return SCS_NULL;
  }

  n_plus_m = A->n + A->m;
  p->m = A->m;
  p->n = A->n;

  p->diag_p = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->diag_r_idxs = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));

  /* Form upper triangular KKT matrix in CSC format */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!p->kkt) {
    scs_printf("Error in forming KKT matrix.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Convert column starts from scs_int (int) to long for Accelerate API */
  p->col_starts = (long *)scs_calloc(n_plus_m + 1, sizeof(long));
  if (!p->col_starts) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  for (i = 0; i <= n_plus_m; i++) {
    p->col_starts[i] = (long)p->kkt->p[i];
  }

  /* Set up Accelerate sparse matrix wrapper pointing to KKT data */
#ifdef SFLOAT
  p->accel_mat = (SparseMatrix_Float){
#else
  p->accel_mat = (SparseMatrix_Double){
#endif
      .structure = {.rowCount = n_plus_m,
                    .columnCount = n_plus_m,
                    .columnStarts = p->col_starts,
                    .rowIndices = p->kkt->i,
                    .attributes =
                        {
                            .kind = SparseSymmetric,
                            .triangle = SparseUpperTriangle,
                        },
                    .blockSize = 1},
      .data = p->kkt->x};

  /* Perform symbolic and numeric factorization */
  p->factorization =
      SparseFactor(SparseFactorizationLDLTUnpivoted, p->accel_mat);
  if (p->factorization.status != SparseStatusOK) {
    scs_printf("Error in Apple Accelerate LDLt factorization: %d\n",
               (int)p->factorization.status);
    SparseCleanup(p->factorization);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  p->factorizations++;

  return p;
}

scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  scs_int n_plus_m = p->n + p->m;
#ifdef SFLOAT
  DenseVector_Float bvec = {.count = n_plus_m, .data = b};
#else
  DenseVector_Double bvec = {.count = n_plus_m, .data = b};
#endif
  SparseSolve(p->factorization, bvec);
  return 0;
}

scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;
  for (i = 0; i < p->n; ++i) {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i) {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }

  /* Refactor with updated values (reuses symbolic factorization) */
  SparseRefactor(p->accel_mat, &(p->factorization));
  if (p->factorization.status != SparseStatusOK) {
    scs_printf("Error in Apple Accelerate LDLt refactorization: %d\n",
               (int)p->factorization.status);
    return -1;
  }
  p->factorizations++;
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    if (p->factorizations > 0) {
      SparseCleanup(p->factorization);
    }
    SCS(cs_spfree)(p->kkt);
    scs_free(p->diag_p);
    scs_free(p->diag_r_idxs);
    scs_free(p->col_starts);
    scs_free(p);
  }
}
