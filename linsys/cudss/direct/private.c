#include "private.h"
#include "linsys.h"

/* In case of error abort freeing p */
#define CUDSS_CHECK_ABORT(call, p, fname)                                      \
  do {                                                                         \
    cudssStatus_t status = call;                                               \
    if (status != CUDSS_STATUS_SUCCESS) {                                      \
      scs_printf("CUDSS call " #fname " returned status = %d\n", status);      \
      scs_free_lin_sys_work(p);                                                \
      return SCS_NULL;                                                         \
    }                                                                          \
  } while (0);

/* In case of error abort freeing p */
#define CUDA_CHECK_ABORT(call, p, fname)                                       \
  do {                                                                         \
    cudaError_t status = call;                                                 \
    if (status != cudaSuccess) {                                               \
      printf("CUDA call " #fname " returned status = %d\n", status);           \
      scs_free_lin_sys_work(p);                                                \
      return SCS_NULL;                                                         \
    }                                                                          \
  } while (0);

/* Return the linear system method name */
const char *scs_get_lin_sys_method() {
  return "sparse-direct-cuDSS";
}

/* Free allocated resources for the linear system solver */
void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    /* Free GPU resources */
    if (p->d_kkt_val)
      cudaFree(p->d_kkt_val);
    if (p->d_kkt_row_ptr)
      cudaFree(p->d_kkt_row_ptr);
    if (p->d_kkt_col_ind)
      cudaFree(p->d_kkt_col_ind);
    if (p->d_b)
      cudaFree(p->d_b);
    if (p->d_sol)
      cudaFree(p->d_sol);

    /* Free cuDSS resources */
    if (p->d_kkt_mat)
      cudssMatrixDestroy(p->d_kkt_mat);
    if (p->d_b_mat)
      cudssMatrixDestroy(p->d_b_mat);
    if (p->d_sol_mat)
      cudssMatrixDestroy(p->d_sol_mat);

    if (p->solver_config)
      cudssConfigDestroy(p->solver_config);
    if (p->solver_data && p->handle)
      cudssDataDestroy(p->handle, p->solver_data);
    if (p->handle)
      cudssDestroy(p->handle);

    /* Free CPU resources */
    if (p->kkt)
      SCS(cs_spfree)(p->kkt);
    if (p->sol)
      scs_free(p->sol);
    if (p->diag_r_idxs)
      scs_free(p->diag_r_idxs);
    if (p->diag_p)
      scs_free(p->diag_p);

    scs_free(p);
  }
}

/* Initialize the linear system solver workspace */
ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  ScsLinSysWork *p = scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p)
    return SCS_NULL;

  /* Store problem dimensions */
  p->n = A->n;
  p->m = A->m;
  p->n_plus_m = p->n + p->m;

  /* Allocate CPU memory */
  p->sol = (scs_float *)scs_malloc(sizeof(scs_float) * p->n_plus_m);
  if (!p->sol) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  p->diag_r_idxs = (scs_int *)scs_calloc(p->n_plus_m, sizeof(scs_int));
  if (!p->diag_r_idxs) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  p->diag_p = (scs_float *)scs_calloc(p->n, sizeof(scs_float));
  if (!p->diag_p) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Form KKT matrix as upper-triangular, CSC */
  /* Because of symmetry it is equivalent to lower-triangular, CSR */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!p->kkt) {
    scs_printf("Error in forming KKT matrix");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  cudssStatus_t status;
  cudaError_t cuda_error;

  /* Create cuDSS handle */
  CUDSS_CHECK_ABORT(cudssCreate(&p->handle), p, "cudssCreate");
  /* Creating cuDSS solver configuration and data objects */

  CUDSS_CHECK_ABORT(cudssConfigCreate(&p->solver_config), p,
                    "cudssConfigCreate");
  CUDSS_CHECK_ABORT(cudssDataCreate(p->handle, &p->solver_data), p,
                    "cudssDataCreate");

  /* Allocate device memory for KKT matrix */
  scs_int nnz = p->kkt->p[p->n_plus_m];

  CUDA_CHECK_ABORT(cudaMalloc((void **)&p->d_kkt_val, nnz * sizeof(scs_float)),
                   p, "cudaMalloc: kkt_val");
  CUDA_CHECK_ABORT(cudaMalloc((void **)&p->d_kkt_row_ptr,
                              (p->n_plus_m + 1) * sizeof(scs_int)),
                   p, "cudaMalloc: kkt_row_ptr");
  CUDA_CHECK_ABORT(
      cudaMalloc((void **)&p->d_kkt_col_ind, nnz * sizeof(scs_int)), p,
      "cudaMalloc: kkt_col_ind");

  /* Copy KKT matrix to device */
  /* Note: we treat column pointers (p->kkt->p) as row pointers on the device */
  CUDA_CHECK_ABORT(cudaMemcpy(p->d_kkt_val, p->kkt->x, nnz * sizeof(scs_float),
                              cudaMemcpyHostToDevice),
                   p, "cudaMemcpy: kkt_val");
  CUDA_CHECK_ABORT(cudaMemcpy(p->d_kkt_row_ptr, p->kkt->p,
                              (p->kkt->n + 1) * sizeof(scs_int),
                              cudaMemcpyHostToDevice),
                   p, "cudaMemcpy: kkt_row_ptr");
  CUDA_CHECK_ABORT(cudaMemcpy(p->d_kkt_col_ind, p->kkt->i,
                              nnz * sizeof(scs_int), cudaMemcpyHostToDevice),
                   p, "cudaMemcpy: kkt_col_ind");

  /* Create kkt matrix descriptor */
  /* We pass the kkt matrix as symmetric, lower triangular */
  cudssMatrixType_t mtype = CUDSS_MTYPE_SYMMETRIC;
  cudssMatrixViewType_t mview = CUDSS_MVIEW_LOWER;
  cudssIndexBase_t base = CUDSS_BASE_ZERO;
  CUDSS_CHECK_ABORT(cudssMatrixCreateCsr(
                        &p->d_kkt_mat, p->kkt->m, p->kkt->n, nnz,
                        p->d_kkt_row_ptr, NULL, p->d_kkt_col_ind, p->d_kkt_val,
                        SCS_CUDA_INDEX, SCS_CUDA_FLOAT, mtype, mview, base),
                    p, "cudssMatrixCreateCsr");

  /* Allocate device memory for vectors */
  CUDA_CHECK_ABORT(
      cudaMalloc((void **)&p->d_b, p->n_plus_m * sizeof(scs_float)), p,
      "cudaMalloc: b");
  CUDA_CHECK_ABORT(
      cudaMalloc((void **)&p->d_sol, p->n_plus_m * sizeof(scs_float)), p,
      "cudaMalloc: sol");

  /* Create RHS and solution matrix descriptors */
  scs_int nrhs = 1;
  CUDSS_CHECK_ABORT(cudssMatrixCreateDn(&p->d_b_mat, p->n_plus_m, nrhs,
                                        p->n_plus_m, p->d_b, SCS_CUDA_FLOAT,
                                        CUDSS_LAYOUT_COL_MAJOR),
                    p, "cudssMatrixCreateDn: b");
  CUDSS_CHECK_ABORT(cudssMatrixCreateDn(&p->d_sol_mat, p->n_plus_m, nrhs,
                                        p->n_plus_m, p->d_sol, SCS_CUDA_FLOAT,
                                        CUDSS_LAYOUT_COL_MAJOR),
                    p, "cudssMatrixCreateDn: sol");

  /* Symbolic factorization */
  CUDSS_CHECK_ABORT(cudssExecute(p->handle, CUDSS_PHASE_ANALYSIS,
                                 p->solver_config, p->solver_data, p->d_kkt_mat,
                                 p->d_sol_mat, p->d_b_mat),
                    p, "cudssExecute: analysis");

  /* Numerical Factorization */
  CUDSS_CHECK_ABORT(cudssExecute(p->handle, CUDSS_PHASE_FACTORIZATION,
                                 p->solver_config, p->solver_data, p->d_kkt_mat,
                                 p->d_sol_mat, p->d_b_mat),
                    p, "cudssExecute: factorization");

  return p;
}

/* Solve the linear system for a given RHS b */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *ws,
                          scs_float tol) {
  /* Copy right-hand side to device */
  cudaError_t custatus = cudaMemcpy(p->d_b, b, p->n_plus_m * sizeof(scs_float),
                                    cudaMemcpyHostToDevice);
  if (custatus != cudaSuccess) {
    scs_printf("scs_solve_lin_sys: Error copying `b` side to device: %d\n",
               (int)custatus);
    return custatus;
  }

  // is this really needed?
  cudssMatrixSetValues(p->d_b_mat, p->d_b);

  /* Solve the system */
  cudssStatus_t status =
      cudssExecute(p->handle, CUDSS_PHASE_SOLVE, p->solver_config,
                   p->solver_data, p->d_kkt_mat, p->d_sol_mat, p->d_b_mat);

  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf("scs_solve_lin_sys: Error during solve: %d\n", (int)status);
    return status;
  }

  /* Copy solution back to host */
  custatus = cudaMemcpy(b, p->d_sol, p->n_plus_m * sizeof(scs_float),
                        cudaMemcpyDeviceToHost);
  if (status != cudaSuccess) {
    scs_printf("scs_solve_lin_sys: Error copying d_sol to host: %d\n",
               (int)status);
    return status;
  }

  return 0; /* Success */
}

/* Update the KKT matrix when R changes */
void scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;

  /* Update KKT matrix on CPU */
  for (i = 0; i < p->n; ++i) {
    /* top left is R_x + P */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i) {
    /* bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }

  /* Copy updated values to device */
  cudaError_t custatus = cudaMemcpy(p->d_kkt_val, p->kkt->x,
                                    p->kkt->p[p->n_plus_m] * sizeof(scs_float),
                                    cudaMemcpyHostToDevice);
  if (custatus != cudaSuccess) {
    scs_printf(
        "scs_update_lin_sys_diag_r: Error copying kkt->x to device: %d\n",
        (int)custatus);
    return;
  }

  /* Update the matrix values in cuDSS */
  cudssStatus_t status;
  status = cudssMatrixSetCsrPointers(p->d_kkt_mat, p->d_kkt_row_ptr, NULL,
                                     p->d_kkt_col_ind, p->d_kkt_val);
  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf(
        "scs_update_lin_sys_diag_r: Error updating kkt matrix on device: %d\n",
        (int)status);
    return;
  }

  /* Perform Refactorization with the updated matrix */
  status =
      cudssExecute(p->handle, CUDSS_PHASE_REFACTORIZATION, p->solver_config,
                   p->solver_data, p->d_kkt_mat, p->d_sol_mat, p->d_b_mat);
  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf("scs_update_lin_sys_diag_r: Error during re-factorization: %d\n",
               (int)status);
    return;
  }
}
