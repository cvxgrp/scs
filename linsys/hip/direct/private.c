#include "private.h"

#include <hip/hip_runtime.h>
#include <hipsparse/hipsparse.h>

const hipsparseOperation_t NON_TRANSP = HIPSPARSE_OPERATION_NON_TRANSPOSE;
const hipsparseSolvePolicy_t NO_LVL_POLICY = HIPSPARSE_SOLVE_POLICY_NO_LEVEL;

const char *scs_get_lin_sys_method()
{
  return "hipsparse-direct";
}

void scs_free_lin_sys_work(ScsLinSysWork *work)
{
  if (work == NULL)
    return;

  // Free device memory
  if (work->d_vals)
    hipFree(work->d_vals);
  if (work->d_row_ptrs)
    hipFree(work->d_row_ptrs);
  if (work->d_col_inds)
    hipFree(work->d_col_inds);
  if (work->d_b)
    hipFree(work->d_b);
  if (work->d_x)
    hipFree(work->d_x);
  if (work->buffer)
    hipFree(work->buffer);

  // Free LU decomposition info
  if (work->info_LU)
    hipsparseDestroyCsrsv2Info(work->info_LU);

  // Free HIPSPARSE matrix descriptor and handle
  if (work->descr)
    hipsparseDestroyMatDescr(work->descr);
  if (work->handle)
    hipsparseDestroy(work->handle);

  // Free the matrix kkt data
  if (work->kkt)
    (SCS(cs_spfree)(work->kkt));

  // Free host-side arrays used for updates
  if (work->diag_r_idxs)
    scs_free(work->diag_r_idxs);
  if (work->diag_p)
    scs_free(work->diag_p);

  // Finally, free the work struct itself
  scs_free(work);
}

hipsparseStatus_t __initialize_work(ScsLinSysWork *work)
{
  ScsMatrix *A = work->kkt;
  hipsparseStatus_t status;

  // Initialize matrix descriptor
  status = hipsparseCreateMatDescr(&(work->descr));
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- descriptor: %d.\n", (int)status);
  }
  hipsparseSetMatIndexBase(work->descr, HIPSPARSE_INDEX_BASE_ZERO); // Zero-based indexing
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- index 0: %d.\n", (int)status);
  }
  status = hipsparseSetMatType(work->descr, HIPSPARSE_MATRIX_TYPE_SYMMETRIC); // Symmetric matrix
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- type symmetric: %d.\n", (int)status);
  }
  status = hipsparseSetMatFillMode(work->descr, HIPSPARSE_FILL_MODE_UPPER); // stored in upper-diagonal part
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- fill upper: %d.\n", (int)status);
  }
  status = hipsparseSetMatDiagType(work->descr, HIPSPARSE_DIAG_TYPE_NON_UNIT); // with non-unit diagonal
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- diagonal non-unit: %d.\n", (int)status);
  }

  // Compute number of non-zeros
  int nnz = A->p[A->n]; // The last element of A->p gives the number of non-zeros

  // Allocate memory on device
  hipMalloc(&(work->d_vals), nnz * sizeof(scs_float));          // Matrix values (number of non-zeros)
  hipMalloc(&(work->d_row_ptrs), (A->n + 1) * sizeof(scs_int)); // Column pointers (size n + 1)
  hipMalloc(&(work->d_col_inds), nnz * sizeof(scs_int));        // Row indices (number of non-zeros)

  // Preallocate memory for vectors b and x (for solving Ax = b)
  hipMalloc(&(work->d_b), A->m * sizeof(scs_float)); // RHS vector b
  hipMalloc(&(work->d_x), A->m * sizeof(scs_float)); // Solution vector x

  // Copy matrix to device
  hipMemcpy(work->d_vals, A->x, nnz * sizeof(scs_float), hipMemcpyHostToDevice);
  hipMemcpy(work->d_row_ptrs, A->p, (A->n + 1) * sizeof(scs_int), hipMemcpyHostToDevice);
  hipMemcpy(work->d_col_inds, A->i, nnz * sizeof(scs_int), hipMemcpyHostToDevice);

  // Initialize HIPSPARSE
  status = hipsparseCreate(&(work->handle));
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- handle: %d.\n", (int)status);
  }
  // Create info object for LU decomposition
  status = hipsparseCreateCsrsv2Info(&(work->info_LU));
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- create info: %d.\n", (int)status);
  }

  // Analyze step to get buffer size
  status = hipsparseDcsrsv2_bufferSize(work->handle, NON_TRANSP, A->m, nnz, work->descr,
                                       work->d_vals, work->d_row_ptrs, work->d_col_inds,
                                       work->info_LU, &(work->bufferSize));

  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- bufferSize: %d.\n", (int)status);
  }

  hipMalloc(&(work->buffer), work->bufferSize);

  // Perform symbolic factorization for LU
  status = hipsparseDcsrsv2_analysis(work->handle, NON_TRANSP, A->m, nnz, work->descr,
                                     work->d_vals, work->d_row_ptrs, work->d_col_inds,
                                     work->info_LU, NO_LVL_POLICY, work->buffer);

  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- analysis: %d.\n", (int)status);
  }

  // Perform numeric factorization (LU decomposition)
  scs_float alpha = 1.0;
  status = hipsparseDcsrsv2_solve(work->handle, NON_TRANSP, A->m, nnz, &alpha, work->descr,
                                  work->d_vals, work->d_row_ptrs, work->d_col_inds,
                                  work->info_LU, NULL, NULL, NO_LVL_POLICY, work->buffer);
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in init -- numeric fact: %d.\n", (int)status);
  }
  return status;
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r)
{
  ScsLinSysWork *p = scs_calloc(1, sizeof(ScsLinSysWork));

  p->n = A->n;
  p->m = A->m;
  scs_int n_plus_m = p->n + p->m;

  p->diag_r_idxs = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->diag_p = (scs_float *)scs_calloc(p->n, sizeof(scs_float));

  // p->kkt is CSC in lower triangular form; this is equivalen to upper CSR
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 0);
  if (!(p->kkt))
  {
    scs_printf("Error in forming KKT matrix");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  hipsparseStatus_t status;
  status = __initialize_work(p);

  if (status == HIPSPARSE_STATUS_SUCCESS)
  {
    return p;
  }
  else
  {
    scs_printf("Error in factorisation: %d.\n", (int)status);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
}

/* Returns solution to linear system Ax = b with solution stored in b */
scs_int scs_solve_lin_sys(ScsLinSysWork *work, scs_float *b, const scs_float *ws,
                          scs_float tol)
{
  if (work == NULL || b == NULL || ws == NULL)
  {
    return -1; // Error: invalid input
  }

  // Copy warmstart solution to device
  hipMemcpy(work->d_x, ws, work->kkt->m * sizeof(scs_float), hipMemcpyHostToDevice);

  // Copy b to device
  hipMemcpy(work->d_b, b, work->kkt->m * sizeof(scs_float), hipMemcpyHostToDevice);

  // Solve the system Ax = b using LU decomposition
  hipsparseStatus_t status;
  scs_float alpha = 1.0;
  scs_int nnz = work->kkt->p[work->kkt->n];
  status = hipsparseDcsrsv2_solve(work->handle, NON_TRANSP, work->kkt->m, nnz, &alpha, work->descr,
                                  work->d_vals, work->d_row_ptrs, work->d_col_inds,
                                  work->info_LU, work->d_x, work->d_b, NO_LVL_POLICY, work->buffer);
  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error during linear system solution: %d.\n", (int)status);
  }

  // Copy the solution back to the host
  hipMemcpy(b, work->d_x, work->kkt->m * sizeof(scs_float), hipMemcpyDeviceToHost);

  return (scs_int)status;
}

/* Update factorization when R changes */
void scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r)
{
  scs_int i;

  for (i = 0; i < p->n; ++i)
  {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i)
  {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }

  scs_int nnz = p->kkt->p[p->kkt->n];
  hipMemcpy(p->d_vals, p->kkt->x, nnz * sizeof(scs_float), hipMemcpyHostToDevice);

  // Perform numeric factorization (LU decomposition) after changes
  hipsparseStatus_t status;
  scs_float alpha = 1.0;
  status = hipsparseDcsrsv2_solve(p->handle, NON_TRANSP, p->kkt->m, nnz, &alpha, p->descr,
                                  p->d_vals, p->d_row_ptrs, p->d_col_inds,
                                  p->info_LU, NULL, NULL, NO_LVL_POLICY, p->buffer);

  if (status != HIPSPARSE_STATUS_SUCCESS)
  {
    scs_printf("Error in factorization when updating: %d.\n",
               (int)status);
    scs_free_lin_sys_work(p);
  }
}
