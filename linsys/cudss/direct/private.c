#include "private.h"
#include "linsys.h"
#include <limits.h>
#include <string.h>

/* The smallest K (nonzeros of the stored triangle, diagonal included) whose
 * products are computed on the GPU; below it scs_lin_sys_residual says "not
 * provided" and SCS computes them on the CPU, from the problem data. A GPU
 * product costs about 0.04 ms whatever the size (the copy of b, kernel
 * launches, the amax and its scalar copy), SCS's own about 1 to 1.5 ns per
 * nonzero. Measured per product on an NVIDIA GB10 (Grace CPU, unified
 * memory, 33 problems): the GPU was 1.3 to 50 times slower up to 4.2e4
 * nonzeros (0.04 against 0.0008 ms at 1.2e3) and faster from 9.0e4 on (0.64
 * to 0.87 times up to 1.2e5, 0.12 to 0.37 times from 3.7e5; BOYD1, 9.3e5:
 * 0.18 against 1.27 ms). In the solves, on random LPs and SOCPs with 1.3e4
 * nonzeros at 1e-9, an iteration with the correction form on cost 1.6 times
 * a plain one with the GPU products against 1.1 with SCS's. A discrete GPU,
 * where copies cross PCIe, likely moves the break-even point up; not
 * measured. Deterministic: the choice depends on the problem's size only.
 * At or above it the products are part of the solver: if they cannot be set
 * up (GPU memory for the second copy of K, a cuBLAS or cuSPARSE error, K too
 * large for the index type, see products_setup), scs_init_lin_sys_work fails
 * with a message, rather than leaving the products to SCS. Setup cost
 * (products_setup, NVIDIA GB10, random K): 3, 35 and 340 ms at 1e5, 1e6
 * and 8e6 nonzeros, 5 to 13% of cuDSS's analysis and factorization; of that
 * cusparseCreate and cublasCreate took 0.4 to 0.8 ms. */
#ifndef CUDSS_MIN_PRODUCT_NNZ
#define CUDSS_MIN_PRODUCT_NNZ (60000)
#endif

static scs_int products_setup(ScsLinSysWork *p);
static scs_int products_refresh(ScsLinSysWork *p);
static void products_free(ScsLinSysWork *p);

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
    /* Free cuDSS resources first, before freeing the GPU memory they reference */
    if (p->solver_data && p->handle)
      cudssDataDestroy(p->handle, p->solver_data);
    if (p->solver_config)
      cudssConfigDestroy(p->solver_config);

    if (p->d_kkt_mat)
      cudssMatrixDestroy(p->d_kkt_mat);
    if (p->d_b_mat)
      cudssMatrixDestroy(p->d_b_mat);
    if (p->d_sol_mat)
      cudssMatrixDestroy(p->d_sol_mat);

    if (p->handle)
      cudssDestroy(p->handle);

    /* Free GPU memory */
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

    /* Free pinned host memory */
    if (p->h_b_pinned)
      cudaFreeHost(p->h_b_pinned);
    if (p->h_sol_pinned)
      cudaFreeHost(p->h_sol_pinned);

    products_free(p);

    /* Free CPU resources */
    if (p->kkt)
      SCS(cs_spfree)(p->kkt);
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

  /* Create cuDSS handle */
  CUDSS_CHECK_ABORT(cudssCreate(&p->handle), p, "cudssCreate");
  /* Creating cuDSS solver configuration and data objects */

  CUDSS_CHECK_ABORT(cudssConfigCreate(&p->solver_config), p,
                    "cudssConfigCreate");
  CUDSS_CHECK_ABORT(cudssDataCreate(p->handle, &p->solver_data), p,
                    "cudssDataCreate");

  /* Plain LDL' with no pivoting, matching or refinement, as in the CPU direct
   * solver (QDLDL): the KKT matrix [R_x + P, A'; A, -R_y] is quasi-definite,
   * so the factorization exists under any ordering. From version 1.0 cuDSS by
   * default applies a matching to a symmetric indefinite matrix. With its
   * default pivoting that makes the factor up to four times larger (2.4 times
   * in geometric mean on 17 LPs and QPs), and without pivoting it made the
   * solves inaccurate enough to stall SCS on some LPs. Tradeoff: one solve
   * per iteration, but the accuracy of a solve rests on the conditioning of
   * the system (worst test-suite residual 6e-8, against 5e-11 with cuDSS's
   * default pivoting and matching).
   *
   * Single precision keeps one refinement step per solve, which about doubles
   * the cost of a solve: without it the test suite fails (rob_gauss_cov_est
   * does not converge), and Bunch-Kaufman pivoting instead of refinement does
   * not fix that. */
  {
    cudssPivotType_t pivot_type = CUDSS_PIVOT_NONE;
#ifdef SFLOAT
    int ir_n_steps = 1;
#else
    int ir_n_steps = 0; /* cuDSS's default, set so that it stays fixed */
#endif
    CUDSS_CHECK_ABORT(cudssConfigSet(p->solver_config, CUDSS_CONFIG_PIVOT_TYPE,
                                     &pivot_type, sizeof(pivot_type)),
                      p, "cudssConfigSet: pivot type");
    CUDSS_CHECK_ABORT(cudssConfigSet(p->solver_config, CUDSS_CONFIG_IR_N_STEPS,
                                     &ir_n_steps, sizeof(ir_n_steps)),
                      p, "cudssConfigSet: refinement steps");
#if defined(CUDSS_VERSION) && CUDSS_VERSION >= 10000
    /* cuDSS 1.0 turned matching on by default (CUDSS_MATCHING_ALG_AUTO).
     * Older versions do not match unless asked to (0.5 has no matching, 0.6
     * and 0.7 default to CUDSS_CONFIG_USE_MATCHING = 0, 0.8 to
     * CUDSS_MATCHING_ALG_NONE), so they need no call. */
    cudssMatchingAlg_t matching_alg = CUDSS_MATCHING_ALG_NONE;
    CUDSS_CHECK_ABORT(cudssConfigSet(p->solver_config,
                                     CUDSS_CONFIG_MATCHING_ALG, &matching_alg,
                                     sizeof(matching_alg)),
                      p, "cudssConfigSet: matching");
#endif
  }

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
  /* cuDSS 0.8.0+ added an offsetType parameter before indexType. The two
   * arrays use the same scs_int type for us, so we pass SCS_CUDA_INDEX
   * twice. */
#if SCS_CUDSS_NEW_API
  CUDSS_CHECK_ABORT(cudssMatrixCreateCsr(
                        &p->d_kkt_mat, p->kkt->m, p->kkt->n, nnz,
                        p->d_kkt_row_ptr, NULL, p->d_kkt_col_ind, p->d_kkt_val,
                        SCS_CUDA_INDEX, SCS_CUDA_INDEX, SCS_CUDA_FLOAT,
                        mtype, mview, base),
                    p, "cudssMatrixCreateCsr");
#else
  CUDSS_CHECK_ABORT(cudssMatrixCreateCsr(
                        &p->d_kkt_mat, p->kkt->m, p->kkt->n, nnz,
                        p->d_kkt_row_ptr, NULL, p->d_kkt_col_ind, p->d_kkt_val,
                        SCS_CUDA_INDEX, SCS_CUDA_FLOAT, mtype, mview, base),
                    p, "cudssMatrixCreateCsr");
#endif

  /* Allocate device memory for vectors */
  CUDA_CHECK_ABORT(
      cudaMalloc((void **)&p->d_b, p->n_plus_m * sizeof(scs_float)), p,
      "cudaMalloc: b");
  CUDA_CHECK_ABORT(
      cudaMalloc((void **)&p->d_sol, p->n_plus_m * sizeof(scs_float)), p,
      "cudaMalloc: sol");

  /* Allocate pinned host memory for faster H<->D transfers in solve loop */
  CUDA_CHECK_ABORT(
      cudaMallocHost((void **)&p->h_b_pinned, p->n_plus_m * sizeof(scs_float)),
      p, "cudaMallocHost: b");
  CUDA_CHECK_ABORT(cudaMallocHost((void **)&p->h_sol_pinned,
                                  p->n_plus_m * sizeof(scs_float)),
                   p, "cudaMallocHost: sol");

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

  /* Inertia check. The KKT matrix has n positive eigenvalues exactly when
   * R_x + P + A' R_y^{-1} A is positive definite, which holds whenever P is
   * PSD. Fewer means P is not PSD (a P that is only mildly indefinite can
   * still pass). cuDSS returns the two counts in the matrix's index type
   * (SCS_CUDA_INDEX): 64-bit integers under DLONG, where it writes 16 bytes,
   * and 32-bit ones otherwise, where it rejects a 16-byte buffer (cuDSS 1.0.0).
   * scs_int follows DLONG in the same way, so it matches; test_non_psd_p
   * covers both widths. */
  {
    scs_int inertia[2]; /* [0] = positive, [1] = negative */
    size_t inertia_written;
    CUDSS_CHECK_ABORT(
        cudssDataGet(p->handle, p->solver_data, CUDSS_DATA_INERTIA, inertia,
                     sizeof(inertia), &inertia_written),
        p, "cudssDataGet: inertia");
    if (inertia[0] < p->n) {
      scs_printf("KKT matrix has < n positive eigenvalues. P not PSD.\n");
      scs_free_lin_sys_work(p);
      return SCS_NULL;
    }
  }

  /* the refinement's products with K on the GPU for a large enough K (see
   * CUDSS_MIN_PRODUCT_NNZ); setup fails if they cannot be set up */
  if (p->kkt->p[p->n_plus_m] >= CUDSS_MIN_PRODUCT_NNZ) {
    if (products_setup(p) != 0) {
      scs_printf("cuDSS: K has %li nonzeros, at least CUDSS_MIN_PRODUCT_NNZ "
                 "= %li, so the refinement's products with K run on the GPU, "
                 "and they could not be set up (see above)\n",
                 (long)p->kkt->p[p->n_plus_m], (long)CUDSS_MIN_PRODUCT_NNZ);
      scs_free_lin_sys_work(p);
      return SCS_NULL;
    }
    p->products_ok = 1;
  }

  return p;
}

/* Solve the linear system for a given RHS b */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *ws,
                          scs_float tol) {
  size_t vec_bytes = p->n_plus_m * sizeof(scs_float);

  /* Copy RHS to pinned staging buffer, then to device */
  memcpy(p->h_b_pinned, b, vec_bytes);
  cudaError_t custatus =
      cudaMemcpy(p->d_b, p->h_b_pinned, vec_bytes, cudaMemcpyHostToDevice);
  if (custatus != cudaSuccess) {
    scs_printf("scs_solve_lin_sys: Error copying b to device: %d\n",
               (int)custatus);
    return custatus;
  }

  /* Solve the system */
  cudssStatus_t status =
      cudssExecute(p->handle, CUDSS_PHASE_SOLVE, p->solver_config,
                   p->solver_data, p->d_kkt_mat, p->d_sol_mat, p->d_b_mat);

  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf("scs_solve_lin_sys: Error during solve: %d\n", (int)status);
    return status;
  }

  /* Copy solution from device to pinned staging buffer, then to output */
  custatus =
      cudaMemcpy(p->h_sol_pinned, p->d_sol, vec_bytes, cudaMemcpyDeviceToHost);
  if (custatus != cudaSuccess) {
    scs_printf("scs_solve_lin_sys: Error copying d_sol to host: %d\n",
               (int)custatus);
    return custatus;
  }
  memcpy(b, p->h_sol_pinned, vec_bytes);

  return 0; /* Success */
}

/* Update the KKT matrix when R changes */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
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

  /* the copy of K with both triangles used for the products */
  if (p->products_ok && products_refresh(p) != 0) {
    scs_printf("scs_update_lin_sys_diag_r: Error updating K for products\n");
    return -1;
  }

  /* Copy updated values to device */
  cudaError_t custatus = cudaMemcpy(p->d_kkt_val, p->kkt->x,
                                    p->kkt->p[p->n_plus_m] * sizeof(scs_float),
                                    cudaMemcpyHostToDevice);
  if (custatus != cudaSuccess) {
    scs_printf(
        "scs_update_lin_sys_diag_r: Error copying kkt->x to device: %d\n",
        (int)custatus);
    return (scs_int)custatus;
  }

  /* Update the matrix values in cuDSS */
  cudssStatus_t status;
  status = cudssMatrixSetCsrPointers(p->d_kkt_mat, p->d_kkt_row_ptr, NULL,
                                     p->d_kkt_col_ind, p->d_kkt_val);
  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf(
        "scs_update_lin_sys_diag_r: Error updating kkt matrix on device: %d\n",
        (int)status);
    return (scs_int)status;
  }

  /* Perform Refactorization with the updated matrix */
  status =
      cudssExecute(p->handle, CUDSS_PHASE_REFACTORIZATION, p->solver_config,
                   p->solver_data, p->d_kkt_mat, p->d_sol_mat, p->d_b_mat);
  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf("scs_update_lin_sys_diag_r: Error during re-factorization: %d\n",
               (int)status);
    return (scs_int)status;
  }
  return 0;
}

/* ------------------------------------------------------------------------
 * Products with K for SCS's per-iteration refinement (optional, linsys.h).
 *
 * r = b - K x is one cuSPARSE SpMV with K stored with both triangles, and
 * ||r||_inf one cuBLAS amax; the guess x0 and r stay on the GPU, so the
 * correction form x = x0 + K^{-1} r costs the copy of b to the GPU and of x
 * back, as a plain solve does, plus the SpMV, amax and an axpy, and no other
 * copies. SpMV uses CUSPARSE_SPMV_CSR_ALG2, which gives the same result
 * every run (CSR_ALG1, which SCS's GPU indirect solver uses, does not), and
 * amax and axpy are reproducible too, so these products add no run-to-run
 * variation; the cuDSS solve itself is not reproducible run to run (see
 * docs, "Reproducibility of the GPU solvers"). Cost: a second copy of K on
 * the GPU (both triangles, about twice the nonzeros of the stored one) and a
 * host copy of its values for updates of R.
 * ---------------------------------------------------------------------- */
#ifndef SFLOAT
#define SCS_CUBLAS_IAMAX cublasIdamax
#define SCS_CUBLAS_AXPY cublasDaxpy
#else
#define SCS_CUBLAS_IAMAX cublasIsamax
#define SCS_CUBLAS_AXPY cublasSaxpy
#endif

static void products_free(ScsLinSysWork *p) {
  if (p->sp_k)
    cusparseDestroySpMat(p->sp_k);
  if (p->sp_r)
    cusparseDestroyDnVec(p->sp_r);
  if (p->sp_x0)
    cusparseDestroyDnVec(p->sp_x0);
  if (p->sp_xv)
    cusparseDestroyDnVec(p->sp_xv);
  if (p->sp_handle)
    cusparseDestroy(p->sp_handle);
  if (p->bl_handle)
    cublasDestroy(p->bl_handle);
  if (p->d_full_val)
    cudaFree(p->d_full_val);
  if (p->d_full_row_ptr)
    cudaFree(p->d_full_row_ptr);
  if (p->d_full_col_ind)
    cudaFree(p->d_full_col_ind);
  if (p->d_r)
    cudaFree(p->d_r);
  if (p->d_x0)
    cudaFree(p->d_x0);
  if (p->d_xv)
    cudaFree(p->d_xv);
  if (p->d_spmv_buf)
    cudaFree(p->d_spmv_buf);
  if (p->full_src)
    scs_free(p->full_src);
  if (p->full_val_h)
    scs_free(p->full_val_h);
  p->sp_k = SCS_NULL;
  p->sp_r = p->sp_x0 = p->sp_xv = SCS_NULL;
  p->sp_handle = SCS_NULL;
  p->bl_handle = SCS_NULL;
  p->d_full_val = p->d_r = p->d_x0 = p->d_xv = SCS_NULL;
  p->d_full_row_ptr = p->d_full_col_ind = SCS_NULL;
  p->d_spmv_buf = SCS_NULL;
  p->full_src = SCS_NULL;
  p->full_val_h = SCS_NULL;
}

/* full_val_h = the values of K (both triangles) from kkt->x, copied to the
 * GPU */
static scs_int products_refresh(ScsLinSysWork *p) {
  scs_int k;
  for (k = 0; k < p->full_nnz; ++k) {
    p->full_val_h[k] = p->kkt->x[p->full_src[k]];
  }
  return cudaMemcpy(p->d_full_val, p->full_val_h,
                    p->full_nnz * sizeof(scs_float),
                    cudaMemcpyHostToDevice) == cudaSuccess
             ? 0
             : -1;
}

/* the largest scs_int, for the size checks of products_setup */
#ifdef DLONG
#define CUDSS_SCS_INT_MAX ((double)LLONG_MAX)
#else
#define CUDSS_SCS_INT_MAX ((double)INT_MAX)
#endif

/* reports why products_setup failed and leaves it */
#define PRODUCTS_FAIL(what)                                                    \
  do {                                                                         \
    scs_printf("cuDSS: setting up the products with K on the GPU: %s\n",       \
               what);                                                          \
    goto done;                                                                 \
  } while (0)

/* Builds K with both triangles from kkt (the lower triangle in CSR: row r
 * holds columns kkt->i[k] <= r), with the columns of each row in increasing
 * order, copies it to the GPU and creates the cuSPARSE and cuBLAS objects
 * (cusparseCreate and cublasCreate are part of this setup cost). Returns 0
 * on success, else prints why and returns -1; the caller frees what was
 * allocated. */
static scs_int products_setup(ScsLinSysWork *p) {
  scs_int nm = p->n_plus_m, r, c, k, q, nnz, ok = -1;
  const scs_int *lp = p->kkt->p, *li = p->kkt->i;
  scs_int *colptr = SCS_NULL, *fill = SCS_NULL, *trow = SCS_NULL;
  scs_int *tsrc = SCS_NULL, *ptr = SCS_NULL, *ind = SCS_NULL;
  size_t vb = nm * sizeof(scs_float), buf_size = 0;
  scs_float alpha = -1., beta = 1.;
  /* sizes: cuBLAS takes the vector length as an int, and K with both
   * triangles (at most twice the stored nonzeros) must fit scs_int, the
   * index type of the cuSPARSE matrix (32-bit unless DLONG) */
  if ((double)nm > (double)INT_MAX) {
    PRODUCTS_FAIL("n + m is larger than cuBLAS's int vector length");
  }
  if (2. * (double)lp[nm] > CUDSS_SCS_INT_MAX) {
    PRODUCTS_FAIL("K with both triangles has too many nonzeros for 32-bit "
                  "indices; build with DLONG=1");
  }
  colptr = (scs_int *)scs_calloc(nm + 1, sizeof(scs_int));
  fill = (scs_int *)scs_calloc(nm + 1, sizeof(scs_int));
  if (!colptr || !fill) {
    PRODUCTS_FAIL("out of host memory");
  }
  /* entries per column of the full matrix (per row too, by symmetry) */
  for (r = 0; r < nm; ++r) {
    for (k = lp[r]; k < lp[r + 1]; ++k) {
      c = li[k];
      colptr[c + 1]++;
      if (c != r) {
        colptr[r + 1]++;
      }
    }
  }
  for (r = 0; r < nm; ++r) {
    colptr[r + 1] += colptr[r];
  }
  nnz = colptr[nm];
  trow = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  tsrc = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  ptr = (scs_int *)scs_calloc(nm + 1, sizeof(scs_int));
  ind = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  p->full_src = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  p->full_val_h = (scs_float *)scs_calloc(nnz, sizeof(scs_float));
  if (!trow || !tsrc || !ptr || !ind || !p->full_src || !p->full_val_h) {
    PRODUCTS_FAIL("out of host memory");
  }
  /* bucket the entries by column (entry (r, c) and its mirror (c, r)) */
  memcpy(fill, colptr, nm * sizeof(scs_int));
  for (r = 0; r < nm; ++r) {
    for (k = lp[r]; k < lp[r + 1]; ++k) {
      c = li[k];
      q = fill[c]++;
      trow[q] = r;
      tsrc[q] = k;
      if (c != r) {
        q = fill[r]++;
        trow[q] = c;
        tsrc[q] = k;
      }
    }
  }
  /* then by row, taking the columns in increasing order: every row's
   * columns come out sorted. The row counts are the column counts. */
  memcpy(ptr, colptr, (nm + 1) * sizeof(scs_int));
  memcpy(fill, colptr, nm * sizeof(scs_int));
  for (c = 0; c < nm; ++c) {
    for (q = colptr[c]; q < colptr[c + 1]; ++q) {
      k = fill[trow[q]]++;
      ind[k] = c;
      p->full_src[k] = tsrc[q];
    }
  }
  p->full_nnz = nnz;
  if (cudaMalloc((void **)&p->d_full_val, nnz * sizeof(scs_float)) !=
          cudaSuccess ||
      cudaMalloc((void **)&p->d_full_row_ptr, (nm + 1) * sizeof(scs_int)) !=
          cudaSuccess ||
      cudaMalloc((void **)&p->d_full_col_ind, nnz * sizeof(scs_int)) !=
          cudaSuccess ||
      cudaMalloc((void **)&p->d_r, vb) != cudaSuccess ||
      cudaMalloc((void **)&p->d_x0, vb) != cudaSuccess ||
      cudaMalloc((void **)&p->d_xv, vb) != cudaSuccess) {
    PRODUCTS_FAIL("out of GPU memory (cudaMalloc failed)");
  }
  if (cudaMemcpy(p->d_full_row_ptr, ptr, (nm + 1) * sizeof(scs_int),
                 cudaMemcpyHostToDevice) != cudaSuccess ||
      cudaMemcpy(p->d_full_col_ind, ind, nnz * sizeof(scs_int),
                 cudaMemcpyHostToDevice) != cudaSuccess ||
      cudaMemset(p->d_x0, 0, vb) != cudaSuccess ||
      products_refresh(p) != 0) {
    PRODUCTS_FAIL("copying K to the GPU failed");
  }
  if (cusparseCreate(&p->sp_handle) != CUSPARSE_STATUS_SUCCESS) {
    PRODUCTS_FAIL("cusparseCreate failed");
  }
  if (cublasCreate(&p->bl_handle) != CUBLAS_STATUS_SUCCESS) {
    PRODUCTS_FAIL("cublasCreate failed");
  }
  if (cusparseCreateCsr(&p->sp_k, nm, nm, nnz, p->d_full_row_ptr,
                        p->d_full_col_ind, p->d_full_val, SCS_CUSPARSE_INDEX,
                        SCS_CUSPARSE_INDEX, CUSPARSE_INDEX_BASE_ZERO,
                        SCS_CUSPARSE_FLOAT) != CUSPARSE_STATUS_SUCCESS ||
      cusparseCreateDnVec(&p->sp_r, nm, p->d_r, SCS_CUSPARSE_FLOAT) !=
          CUSPARSE_STATUS_SUCCESS ||
      cusparseCreateDnVec(&p->sp_x0, nm, p->d_x0, SCS_CUSPARSE_FLOAT) !=
          CUSPARSE_STATUS_SUCCESS ||
      cusparseCreateDnVec(&p->sp_xv, nm, p->d_xv, SCS_CUSPARSE_FLOAT) !=
          CUSPARSE_STATUS_SUCCESS) {
    PRODUCTS_FAIL("creating the cuSPARSE descriptors failed");
  }
  if (cusparseSpMV_bufferSize(p->sp_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                              &alpha, p->sp_k, p->sp_xv, &beta, p->sp_r,
                              SCS_CUSPARSE_FLOAT, CUSPARSE_SPMV_CSR_ALG2,
                              &buf_size) != CUSPARSE_STATUS_SUCCESS) {
    PRODUCTS_FAIL("cusparseSpMV_bufferSize failed");
  }
  if (buf_size > 0 && cudaMalloc(&p->d_spmv_buf, buf_size) != cudaSuccess) {
    PRODUCTS_FAIL("out of GPU memory (cudaMalloc failed)");
  }
  ok = 0;
done:
  scs_free(colptr);
  scs_free(fill);
  scs_free(trow);
  scs_free(tsrc);
  scs_free(ptr);
  scs_free(ind);
  return ok;
}

/* r = b - K x on the GPU (x = the guess x0 when x is SCS_NULL); r stays on
 * the GPU in d_r. Returns ||r||_inf, SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED below
 * CUDSS_MIN_PRODUCT_NNZ, or SCS_LIN_SYS_PRODUCTS_FAILED on error. Uses d_xv,
 * d_r and the pinned staging buffers, never d_sol, so the latest solve's
 * solution stays for scs_lin_sys_keep. */
scs_float scs_lin_sys_residual(ScsLinSysWork *p, const scs_float *b,
                               const scs_float *x) {
  size_t vb = p->n_plus_m * sizeof(scs_float);
  scs_float alpha = -1., beta = 1., res = 0.;
  int imax = 0;
  if (!p->products_ok) {
    return SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED;
  }
  /* d_r = b, through the pinned staging buffer of the solve */
  memcpy(p->h_b_pinned, b, vb);
  if (cudaMemcpy(p->d_r, p->h_b_pinned, vb, cudaMemcpyHostToDevice) !=
      cudaSuccess) {
    return SCS_LIN_SYS_PRODUCTS_FAILED;
  }
  if (x) {
    memcpy(p->h_sol_pinned, x, vb);
    if (cudaMemcpy(p->d_xv, p->h_sol_pinned, vb, cudaMemcpyHostToDevice) !=
        cudaSuccess) {
      return SCS_LIN_SYS_PRODUCTS_FAILED;
    }
  }
  /* d_r = -K x + d_r, then the index of its largest |entry| (1-based) */
  if (cusparseSpMV(p->sp_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                   p->sp_k, x ? p->sp_xv : p->sp_x0, &beta, p->sp_r,
                   SCS_CUSPARSE_FLOAT, CUSPARSE_SPMV_CSR_ALG2,
                   p->d_spmv_buf) != CUSPARSE_STATUS_SUCCESS ||
      SCS_CUBLAS_IAMAX(p->bl_handle, (int)p->n_plus_m, p->d_r, 1, &imax) !=
          CUBLAS_STATUS_SUCCESS) {
    return SCS_LIN_SYS_PRODUCTS_FAILED;
  }
  if (imax > 0 && cudaMemcpy(&res, p->d_r + (imax - 1), sizeof(scs_float),
                             cudaMemcpyDeviceToHost) != cudaSuccess) {
    return SCS_LIN_SYS_PRODUCTS_FAILED;
  }
  return ABS(res);
}

/* x0 = the solution of the latest solve, which is still in d_sol (x is that
 * solution, so it is not read) */
void scs_lin_sys_keep(ScsLinSysWork *p, const scs_float *x) {
  (void)x;
  if (p->products_ok) {
    cudaMemcpy(p->d_x0, p->d_sol, p->n_plus_m * sizeof(scs_float),
               cudaMemcpyDeviceToDevice);
  }
}

/* x = x0 + K^{-1} r with r = d_r from the latest scs_lin_sys_residual, all on
 * the GPU; one copy of x to the host */
scs_int scs_lin_sys_correct(ScsLinSysWork *p, scs_float *x) {
  size_t vb = p->n_plus_m * sizeof(scs_float);
  scs_float one = 1.;
  cudssStatus_t status;
  if (!p->products_ok) {
    return -1;
  }
  if (cudaMemcpy(p->d_b, p->d_r, vb, cudaMemcpyDeviceToDevice) !=
      cudaSuccess) {
    return -1;
  }
  status = cudssExecute(p->handle, CUDSS_PHASE_SOLVE, p->solver_config,
                        p->solver_data, p->d_kkt_mat, p->d_sol_mat,
                        p->d_b_mat);
  if (status != CUDSS_STATUS_SUCCESS) {
    scs_printf("scs_lin_sys_correct: Error during solve: %d\n", (int)status);
    return status;
  }
  /* d_sol = x0 + d: one rounding, the same value as SCS's x0 + d on the
   * host */
  if (SCS_CUBLAS_AXPY(p->bl_handle, (int)p->n_plus_m, &one, p->d_x0, 1,
                      p->d_sol, 1) != CUBLAS_STATUS_SUCCESS) {
    return -1;
  }
  if (cudaMemcpy(p->h_sol_pinned, p->d_sol, vb, cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    return -1;
  }
  memcpy(x, p->h_sol_pinned, vb);
  return 0;
}
