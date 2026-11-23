#include "private.h"
#include "linsys.h"

#include <string.h>

/* norm to use when deciding convergence */
/* should be consistent with CG_NORM in glbopts.h */
#define USE_L2_NORM (0)

static scs_float cg_gpu_norm(cublasHandle_t cublas_handle, const scs_float *r,
                             scs_int n) {
#if USE_L2_NORM > 0
  scs_float nrm;
  CUBLAS(nrm2)(cublas_handle, n, r, 1, &nrm);
#else
  int idx; /* cuBLAS's index type, also with DLONG */
  scs_float nrm;
  CUBLASI(amax)(cublas_handle, n, r, 1, &idx);
  /* NOTE: we take idx -1 here since the routine above returns Fortran idxs */
  cudaMemcpy(&nrm, &(r[idx - 1]), sizeof(scs_float), cudaMemcpyDeviceToHost);
  nrm = ABS(nrm);
#endif
  return nrm;
}

const char *scs_get_lin_sys_method() {
  return "sparse-indirect GPU";
}

/* Not possible to do this on the fly due to M_ii += a_i' (R_y)^-1 a_i */
/* set M = inv ( diag ( R_x + P + A' R_y^{-1} A ) ), and with deflation on the
 * diagonal itself in Md_gpu (the filter of a harvest). Returns 0, or -1 if a
 * copy to the GPU failed. */
static scs_int set_preconditioner(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i, k;
  const ScsMatrix *A = p->A;
  const ScsMatrix *P = p->P;
  scs_float *M = p->M;

#if VERBOSITY > 0
  scs_printf("getting pre-conditioner\n");
#endif

  /* M_ii = (R_x)_i + P_ii + a_i' (R_y)^-1 a_i */
  for (i = 0; i < A->n; ++i) { /* cols */
    /* M_ii = (R_x)_i */
    M[i] = diag_r[i];
    /* M_ii += a_i' (R_y)^-1 a_i */
    for (k = A->p[i]; k < A->p[i + 1]; ++k) {
      /* A->i[k] is row of entry k with value A->x[k] */
      M[i] += A->x[k] * A->x[k] / diag_r[A->n + A->i[k]];
    }
    if (P) {
      for (k = P->p[i]; k < P->p[i + 1]; k++) {
        /* diagonal element only */
        if (P->i[k] == i) { /* row == col */
          /* M_ii += P_ii */
          M[i] += P->x[k];
          break;
        }
      }
    }
  }
  if (p->Md_gpu && cudaMemcpy(p->Md_gpu, M, A->n * sizeof(scs_float),
                              cudaMemcpyHostToDevice) != cudaSuccess) {
    return -1;
  }
  /* finally invert for pre-conditioner */
  for (i = 0; i < A->n; ++i) {
    M[i] = 1. / M[i];
  }
  if (cudaMemcpy(p->M_gpu, M, A->n * sizeof(scs_float),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    return -1;
  }
#if VERBOSITY > 0
  scs_printf("finished getting pre-conditioner\n");
#endif
  return 0;
}

scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;

  /* R_x to gpu */
  if (cudaMemcpy(p->r_x_gpu, diag_r, p->n * sizeof(scs_float),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    scs_printf("ERROR: GPU indirect solver: copying R to the GPU failed\n");
    return -1;
  }

  /* 1/R_y to gpu */
  for (i = 0; i < p->m; ++i)
    p->inv_r_y[i] = 1. / diag_r[p->n + i];
  if (cudaMemcpy(p->inv_r_y_gpu, p->inv_r_y, p->m * sizeof(scs_float),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    scs_printf("ERROR: GPU indirect solver: copying R to the GPU failed\n");
    return -1;
  }

  /* set preconditioner M on gpu */
  if (set_preconditioner(p, diag_r) != 0) {
    scs_printf("ERROR: GPU indirect solver: copying the preconditioner to "
               "the GPU failed\n");
    return -1;
  }
  /* the operator and preconditioner changed: the harvested eigenvectors
   * live in the old preconditioned geometry, discard them; the next deep
   * solve (triggered by this same metric change) re-harvests */
  SCS(eigcg_reset)(&p->dfl);
  return 0;
}

/* safe on a partly set up p (scs_init_lin_sys_work calls it on failure) */
void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    SCS(eigcg_free)(&p->dfl);
    scs_free(p->M);
    scs_free(p->inv_r_y);
    cudaFree(p->p);
    cudaFree(p->r);
    cudaFree(p->Gp);
    cudaFree(p->bg);
    cudaFree(p->tmp_m);
    cudaFree(p->z);
    cudaFree(p->M_gpu);
    cudaFree(p->Md_gpu);
    cudaFree(p->dfl_small);
    cudaFree(p->r_x_gpu);
    cudaFree(p->inv_r_y_gpu);
    if (p->Pg) {
      SCS(free_gpu_matrix)(p->Pg);
      scs_free(p->Pg);
    }
    if (p->Ag) {
      SCS(free_gpu_matrix)(p->Ag);
      scs_free(p->Ag);
    }
    if (p->Agt) {
      SCS(free_gpu_matrix)(p->Agt);
      scs_free(p->Agt);
    }
    if (p->buffer != SCS_NULL) {
      cudaFree(p->buffer);
    }
    if (p->dn_vec_m) {
      cusparseDestroyDnVec(p->dn_vec_m);
    }
    if (p->dn_vec_n) {
      cusparseDestroyDnVec(p->dn_vec_n);
    }
    if (p->dn_vec_n_p) {
      cusparseDestroyDnVec(p->dn_vec_n_p);
    }
    if (p->cusparse_handle) {
      cusparseDestroy(p->cusparse_handle);
    }
    if (p->cublas_handle) {
      cublasDestroy(p->cublas_handle);
    }
    /* Don't reset because it interferes with other GPU programs. */
    /* cudaDeviceReset(); */
    SCS(lin_sys_shared_free)(&p->sh);
    scs_free(p);
  }
}

/* z = M * z elementwise in place, assumes M, z on GPU.
 *
 * dgmm computes C = diag(M) * A; here A and C are both z, as an n x 1 matrix.
 * cuBLAS allows C == A when their leading dimensions are equal.
 *
 * This runs three times per CG iteration. It used to be tbmv with zero
 * bandwidth, which gives the same values but whose time grows as n^2 (cuBLAS
 * 13.1: 0.7 ms per call at n = 1e5 and 63 ms at n = 1e6, against 0.01 ms and
 * 0.1 ms for dgmm).
 */
static void scale_by_diag(cublasHandle_t cublas_handle, scs_float *M,
                          scs_float *z, scs_int n) {
  if (n <= 0) {
    return; /* dgmm rejects a leading dimension of 0 */
  }
  CUBLAS(dgmm)
  (cublas_handle, CUBLAS_SIDE_LEFT, n, 1, z, n, M, 1, z, n);
}

/* y = (R_x + P + A' R_y^{-1} A) x */
static void mat_vec(ScsLinSysWork *p, const scs_float *x, scs_float *y) {
  /* x and y MUST already be loaded to GPU */
  scs_float *z = p->tmp_m; /* temp memory */
  cudaMemset(z, 0, p->m * sizeof(scs_float));

  cusparseDnVecSetValues(p->dn_vec_m, (void *)z);
  cusparseDnVecSetValues(p->dn_vec_n, (void *)x);
  cusparseDnVecSetValues(p->dn_vec_n_p, (void *)y);

  /* y = x */
  cudaMemcpy(y, x, p->n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
  /* y = R_x * x */
  scale_by_diag(p->cublas_handle, p->r_x_gpu, y, p->n);

  if (p->Pg) {
    /* y = R_x * x + P x */
    SCS(accum_by_p_gpu)(p->Pg, p->dn_vec_n, p->dn_vec_n_p, p->cusparse_handle,
                        &p->buffer_size, &p->buffer);
  }

  /* z = Ax */
#if GPU_TRANSPOSE_MAT > 0
  SCS(accum_by_atrans_gpu)(p->Agt, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle,
                           &p->buffer_size, &p->buffer);
#else
  SCS(accum_by_a_gpu)(p->Ag, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle,
                      &p->buffer_size, &p->buffer);
#endif
  /* z = R_y^{-1} A x */
  scale_by_diag(p->cublas_handle, p->inv_r_y_gpu, z, p->m);

  /* y += A'z => y = R_x * x + P x + A' R_y^{-1} Ax */
  SCS(accum_by_atrans_gpu)(p->Ag, p->dn_vec_m, p->dn_vec_n_p,
                           p->cusparse_handle, &p->buffer_size, &p->buffer);
}

/* P comes in upper triangular, expand to full
 * First compute triplet version of full matrix, then compress to CSC
 * Returns SCS_NULL if out of memory.
 * */
static ScsMatrix *fill_p_matrix(const ScsMatrix *P) {
  scs_int i, j, k, kk;
  scs_int Pnzmax = 2 * P->p[P->n]; /* upper bound */
  ScsMatrix *P_tmp = SCS(cs_spalloc)(P->n, P->n, Pnzmax, 1, 1);
  ScsMatrix *P_full;
  if (!P_tmp) {
    return SCS_NULL;
  }
  kk = 0;
  for (j = 0; j < P->n; j++) { /* cols */
    for (k = P->p[j]; k < P->p[j + 1]; k++) {
      i = P->i[k]; /* row */
      if (i > j) { /* only upper triangular needed */
        break;
      }
      P_tmp->i[kk] = i;
      P_tmp->p[kk] = j;
      P_tmp->x[kk] = P->x[k];
      kk++;
      if (i == j) { /* diagonal */
        continue;
      }
      P_tmp->i[kk] = j;
      P_tmp->p[kk] = i;
      P_tmp->x[kk] = P->x[k];
      kk++;
    }
  }
  P_full = SCS(cs_compress)(P_tmp, kk, SCS_NULL);
  SCS(cs_spfree)(P_tmp);
  return P_full;
}

/* *ptr = count numbers of size bytes each in GPU memory; returns 0, or -1
 * with a message naming what could not be allocated */
static scs_int gpu_malloc(void **ptr, size_t count, size_t size,
                          const char *what) {
  cudaError_t err = cudaMalloc(ptr, count * size);
  if (err != cudaSuccess) {
    *ptr = SCS_NULL;
    scs_printf("ERROR: GPU indirect solver: could not allocate %s on the GPU "
               "(%.3g MB): %s\n",
               what, (double)(count * size) / 1e6, cudaGetErrorString(err));
    /* setup fails; clear the error so that no later check of
     * cudaGetLastError, in this process, takes it for a failure of other
     * work */
    cudaGetLastError();
    return -1;
  }
  return 0;
}

/* ======================== eigCG deflation ======================== */

/* The n-sized operations of the shared eigCG (linsys/eigcg.h), on the
 * GPU: W, AW and the window are in GPU memory. The small matrices and
 * coefficients of the shared code, on the host, go through p->dfl_small.
 * cuBLAS runs in its default host pointer mode, as in CG, so a dot product
 * or a copy to the host waits for the work before it. The status of every
 * cuBLAS call and copy of these ops is checked, and a failure fails the
 * solve (linsys/eigcg.h). Two parts are CG's own and unchecked, as they are
 * in CG: the product by G in gpu_eig_probe (mat_vec, which checks none of
 * its cuSPARSE calls, copies or scaling) and the norm (cg_gpu_norm, which
 * checks neither its amax nor its copy). */

/* called by SCS(eigcg_init) only, at setup */
static scs_float *gpu_eig_alloc(ScsLinSysWork *p, scs_int rows, scs_int cols) {
  void *v = SCS_NULL;
  size_t count = (size_t)rows * (size_t)cols;
  (void)p;
  if (gpu_malloc(&v, count, sizeof(scs_float), "eigCG deflation's vectors")) {
    return SCS_NULL;
  }
  if (cudaMemset(v, 0, count * sizeof(scs_float)) != cudaSuccess) {
    scs_printf("ERROR: GPU indirect solver: could not clear eigCG deflation's "
               "vectors on the GPU\n");
    cudaFree(v);
    cudaGetLastError();
    return SCS_NULL;
  }
  return (scs_float *)v;
}

static void gpu_eig_free(ScsLinSysWork *p, scs_float *v) {
  (void)p;
  cudaFree(v);
}

static scs_int gpu_eig_scale_copy(ScsLinSysWork *p, scs_float *dst,
                                  const scs_float *src, scs_float a) {
  int n = (int)p->n;
  if (CUBLAS_FN(copy)(p->cublas_handle, n, src, 1, dst, 1) !=
      CUBLAS_STATUS_SUCCESS) {
    return -1;
  }
  return CUBLAS_FN(scal)(p->cublas_handle, n, &a, dst, 1) !=
         CUBLAS_STATUS_SUCCESS;
}

static scs_int gpu_eig_copy(ScsLinSysWork *p, scs_float *dst,
                            const scs_float *src, scs_int cols) {
  return cudaMemcpy(dst, src, (size_t)cols * p->n * sizeof(scs_float),
                    cudaMemcpyDeviceToDevice) != cudaSuccess;
}

/* room for count numbers in p->dfl_small, or SCS_NULL (never: its size,
 * set at init, covers what the shared code asks) */
static scs_float *gpu_eig_small(ScsLinSysWork *p, size_t count) {
  return count <= p->dfl_small_n ? p->dfl_small : SCS_NULL;
}

static scs_int gpu_eig_mul_small(ScsLinSysWork *p, scs_float *C,
                                 const scs_float *A, scs_int k,
                                 const scs_float *B, scs_int cols) {
  scs_float one = 1., zero = 0.;
  int n = (int)p->n;
  size_t nb = (size_t)k * cols;
  scs_float *Bd = gpu_eig_small(p, nb);
  if (!Bd || cudaMemcpy(Bd, B, nb * sizeof(scs_float),
                        cudaMemcpyHostToDevice) != cudaSuccess) {
    return -1;
  }
  return CUBLAS_FN(gemm)(p->cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, n,
                         (int)cols, (int)k, &one, A, n, Bd, (int)k, &zero, C,
                         n) != CUBLAS_STATUS_SUCCESS;
}

/* the filter's measure, as the CPU solver's: D v by dgmm with the diagonal
 * (Md_gpu), where the CPU divides by M; they differ at rounding level */
static scs_int gpu_eig_probe(ScsLinSysWork *p, const scs_float *v,
                             scs_float *Kv, scs_float theta, scs_float *rr,
                             scs_float *mm) {
  int n = (int)p->n;
  scs_float neg_onef = -1.;
  scs_float *t = p->z, *d = p->Gp; /* CG's scratch, free after CG */
  cublasHandle_t h = p->cublas_handle;
  mat_vec(p, v, Kv);
  /* t = theta D v, d = Kv - t */
  if (cudaMemcpy(t, v, n * sizeof(scs_float), cudaMemcpyDeviceToDevice) !=
          cudaSuccess ||
      CUBLAS_FN(dgmm)(h, CUBLAS_SIDE_LEFT, n, 1, t, n, p->Md_gpu, 1, t, n) !=
          CUBLAS_STATUS_SUCCESS ||
      CUBLAS_FN(scal)(h, n, &theta, t, 1) != CUBLAS_STATUS_SUCCESS ||
      CUBLAS_FN(dot)(h, n, t, 1, t, 1, mm) != CUBLAS_STATUS_SUCCESS ||
      cudaMemcpy(d, Kv, n * sizeof(scs_float), cudaMemcpyDeviceToDevice) !=
          cudaSuccess ||
      CUBLAS_FN(axpy)(h, n, &neg_onef, t, 1, d, 1) != CUBLAS_STATUS_SUCCESS ||
      CUBLAS_FN(dot)(h, n, d, 1, d, 1, rr) != CUBLAS_STATUS_SUCCESS) {
    return -1;
  }
  return 0;
}

static scs_int gpu_eig_gram(ScsLinSysWork *p, const scs_float *W,
                            const scs_float *AW, scs_int k, scs_float *G) {
  scs_float one = 1., zero = 0.;
  int n = (int)p->n;
  size_t nk = (size_t)k * k;
  scs_float *Gd = gpu_eig_small(p, nk);
  if (!Gd ||
      CUBLAS_FN(gemm)(p->cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, (int)k,
                      (int)k, n, &one, W, n, AW, n, &zero, Gd,
                      (int)k) != CUBLAS_STATUS_SUCCESS) {
    return -1;
  }
  return cudaMemcpy(G, Gd, nk * sizeof(scs_float), cudaMemcpyDeviceToHost) !=
         cudaSuccess;
}

static scs_int gpu_eig_coeffs(ScsLinSysWork *p, const scs_float *W, scs_int k,
                              const scs_float *r, scs_float *c) {
  scs_float one = 1., zero = 0.;
  int n = (int)p->n;
  scs_float *cd = gpu_eig_small(p, (size_t)k);
  if (!cd ||
      CUBLAS_FN(gemv)(p->cublas_handle, CUBLAS_OP_T, n, (int)k, &one, W, n, r,
                      1, &zero, cd, 1) != CUBLAS_STATUS_SUCCESS) {
    return -1;
  }
  return cudaMemcpy(c, cd, k * sizeof(scs_float), cudaMemcpyDeviceToHost) !=
         cudaSuccess;
}

static scs_int gpu_eig_project(ScsLinSysWork *p, const scs_float *W,
                               const scs_float *AW, scs_int k,
                               const scs_float *c, scs_float *x, scs_float *r,
                               scs_int undo) {
  scs_float one = 1., sx = undo ? -1. : 1., sr = undo ? 1. : -1.;
  int n = (int)p->n;
  scs_float *cd = gpu_eig_small(p, (size_t)k);
  if (!cd || cudaMemcpy(cd, c, k * sizeof(scs_float),
                        cudaMemcpyHostToDevice) != cudaSuccess) {
    return -1;
  }
  if (CUBLAS_FN(gemv)(p->cublas_handle, CUBLAS_OP_N, n, (int)k, &sx, W, n, cd,
                      1, &one, x, 1) != CUBLAS_STATUS_SUCCESS) {
    return -1;
  }
  return CUBLAS_FN(gemv)(p->cublas_handle, CUBLAS_OP_N, n, (int)k, &sr, AW, n,
                         cd, 1, &one, r, 1) != CUBLAS_STATUS_SUCCESS;
}

/* CG's norm. A residual with NaN entries: whether cuBLAS's amax picks one
 * (and so whether the projection is undone) is not specified, as with a
 * BLAS's idamax on the CPU; the solve for g's finiteness check catches NaN
 * either way. */
static scs_float gpu_eig_norm(ScsLinSysWork *p, const scs_float *r) {
  return cg_gpu_norm(p->cublas_handle, r, p->n);
}

static const ScsEigcgOps gpu_eig_ops = {
    gpu_eig_alloc,     gpu_eig_free,  gpu_eig_scale_copy, gpu_eig_copy,
    gpu_eig_mul_small, gpu_eig_probe, gpu_eig_gram,       gpu_eig_coeffs,
    gpu_eig_project,   gpu_eig_norm};

/* ======================== Setup ======================== */

static scs_int gpu_upload(void *dst, const void *src, size_t bytes,
                          const char *what) {
  cudaError_t err = cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    scs_printf("ERROR: GPU indirect solver: could not copy %s to the GPU: %s\n",
               what, cudaGetErrorString(err));
    return -1;
  }
  return 0;
}

/* the CSC matrix M (m x n) on the GPU, which cuSPARSE reads as the CSR of
 * M' (n x m); SCS_NULL, with a message, on failure */
static ScsGpuMatrix *gpu_matrix(const ScsMatrix *M, const char *what) {
  scs_int nnz = M->p[M->n];
  ScsGpuMatrix *g = (ScsGpuMatrix *)scs_calloc(1, sizeof(ScsGpuMatrix));
  if (!g) {
    scs_printf("ERROR: GPU indirect solver: out of host memory\n");
    return SCS_NULL;
  }
  g->n = M->n;
  g->m = M->m;
  g->nnz = nnz;
  g->descr = 0;
  if (gpu_malloc((void **)&g->i, nnz, sizeof(scs_int), what) ||
      gpu_malloc((void **)&g->p, M->n + 1, sizeof(scs_int), what) ||
      gpu_malloc((void **)&g->x, nnz, sizeof(scs_float), what) ||
      gpu_upload(g->i, M->i, nnz * sizeof(scs_int), what) ||
      gpu_upload(g->p, M->p, (M->n + 1) * sizeof(scs_int), what) ||
      gpu_upload(g->x, M->x, nnz * sizeof(scs_float), what)) {
    SCS(free_gpu_matrix)(g);
    scs_free(g);
    return SCS_NULL;
  }
  if (cusparseCreateCsr(&g->descr, g->n, g->m, g->nnz, g->p, g->i, g->x,
                        SCS_CUSPARSE_INDEX, SCS_CUSPARSE_INDEX,
                        CUSPARSE_INDEX_BASE_ZERO,
                        SCS_CUDA_FLOAT) != CUSPARSE_STATUS_SUCCESS) {
    scs_printf("ERROR: GPU indirect solver: cuSPARSE could not describe %s\n",
               what);
    g->descr = 0;
    SCS(free_gpu_matrix)(g);
    scs_free(g);
    return SCS_NULL;
  }
  return g;
}

/* grows cuSPARSE's work buffer to at least size bytes */
static scs_int reserve_buffer(ScsLinSysWork *p, size_t size) {
  if (size > p->buffer_size) {
    cudaFree(p->buffer);
    p->buffer = SCS_NULL;
    p->buffer_size = 0;
    if (gpu_malloc(&p->buffer, size, 1, "cuSPARSE's work buffer")) {
      return -1;
    }
    p->buffer_size = size;
  }
  return 0;
}

/* the size of cuSPARSE's work buffer for y += op(A) x, in *size; returns 0
 * or -1 */
static scs_int spmv_buffer_size(ScsLinSysWork *p, cusparseOperation_t op,
                                const ScsGpuMatrix *A, cusparseDnVecDescr_t x,
                                cusparseDnVecDescr_t y, size_t *size) {
  const scs_float onef = 1.0;
  *size = 0;
  return cusparseSpMV_bufferSize(p->cusparse_handle, op, &onef, A->descr, x,
                                 &onef, y, SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
                                 size) != CUSPARSE_STATUS_SUCCESS;
}

/* sizes cuSPARSE's work buffer for every product of the solver (A'z, A x and
 * P x) now, so that the products never allocate (accum_by_*_gpu only grow it
 * when it is too small) and deflation's memory, allocated after this, cannot
 * take the products' */
static scs_int reserve_spmv_buffer(ScsLinSysWork *p) {
  size_t s_at = 0, s_a = 0, s_p = 0;
  scs_int fail = spmv_buffer_size(p, CUSPARSE_OPERATION_NON_TRANSPOSE, p->Ag,
                                  p->dn_vec_m, p->dn_vec_n, &s_at);
#if GPU_TRANSPOSE_MAT > 0
  fail = fail || spmv_buffer_size(p, CUSPARSE_OPERATION_NON_TRANSPOSE, p->Agt,
                                  p->dn_vec_n, p->dn_vec_m, &s_a);
#else
  fail = fail || spmv_buffer_size(p, CUSPARSE_OPERATION_TRANSPOSE, p->Ag,
                                  p->dn_vec_n, p->dn_vec_m, &s_a);
#endif
  if (p->Pg) {
    fail = fail || spmv_buffer_size(p, CUSPARSE_OPERATION_NON_TRANSPOSE, p->Pg,
                                    p->dn_vec_n, p->dn_vec_n_p, &s_p);
  }
  if (fail) {
    scs_printf("ERROR: GPU indirect solver: cuSPARSE could not size its work "
               "buffer\n");
    return -1;
  }
  return reserve_buffer(p, MAX(s_at, MAX(s_a, s_p)));
}

#if GPU_TRANSPOSE_MAT > 0
/* Agt: the transpose of Ag on the GPU, for faster multiplies by A */
static scs_int transpose_a(ScsLinSysWork *p) {
  const ScsMatrix *A = p->A;
  ScsGpuMatrix *Ag = p->Ag, *Agt;
  size_t new_buffer_size = 0;
  scs_int nnz = A->p[A->n];
  Agt = p->Agt = (ScsGpuMatrix *)scs_calloc(1, sizeof(ScsGpuMatrix));
  if (!Agt) {
    scs_printf("ERROR: GPU indirect solver: out of host memory\n");
    return -1;
  }
  Agt->n = A->m;
  Agt->m = A->n;
  Agt->nnz = nnz;
  Agt->descr = 0;
  /* Matrix description */
  if (gpu_malloc((void **)&Agt->i, nnz, sizeof(scs_int), "A's transpose") ||
      gpu_malloc((void **)&Agt->p, A->m + 1, sizeof(scs_int),
                 "A's transpose") ||
      gpu_malloc((void **)&Agt->x, nnz, sizeof(scs_float), "A's transpose")) {
    return -1;
  }
  /* transpose Ag into Agt for faster multiplies */
  /* TODO: memory intensive, could perform transpose in CPU and copy to GPU */
  if (cusparseCsr2cscEx2_bufferSize(
          p->cusparse_handle, A->n, A->m, nnz, Ag->x, Ag->p, Ag->i, Agt->x,
          Agt->p, Agt->i, SCS_CUDA_FLOAT, CUSPARSE_ACTION_NUMERIC,
          CUSPARSE_INDEX_BASE_ZERO, SCS_CSR2CSC_ALG,
          &new_buffer_size) != CUSPARSE_STATUS_SUCCESS ||
      reserve_buffer(p, new_buffer_size) ||
      cusparseCsr2cscEx2(p->cusparse_handle, A->n, A->m, nnz, Ag->x, Ag->p,
                         Ag->i, Agt->x, Agt->p, Agt->i, SCS_CUDA_FLOAT,
                         CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO,
                         SCS_CSR2CSC_ALG,
                         p->buffer) != CUSPARSE_STATUS_SUCCESS ||
      cusparseCreateCsr(&Agt->descr, Agt->n, Agt->m, Agt->nnz, Agt->p, Agt->i,
                        Agt->x, SCS_CUSPARSE_INDEX, SCS_CUSPARSE_INDEX,
                        CUSPARSE_INDEX_BASE_ZERO,
                        SCS_CUDA_FLOAT) != CUSPARSE_STATUS_SUCCESS) {
    scs_printf("ERROR: GPU indirect solver: cuSPARSE could not transpose A\n");
    Agt->descr = 0; /* not created */
    return -1;
  }
  return 0;
}
#endif

static scs_int plain_solve(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol);

/* everything scs_init_lin_sys_work sets up; returns 0, or -1 with a message
 * (the caller then frees p) */
static scs_int init_work(ScsLinSysWork *p, const scs_float *diag_r) {
  cudaError_t err;
  const ScsMatrix *A = p->A, *P = p->P;
  scs_int n = p->n, m = p->m;
  const char *vecs = "the conjugate gradient vectors";

  p->inv_r_y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  p->M = (scs_float *)scs_calloc(n, sizeof(scs_float));
  if (!p->inv_r_y || !p->M) {
    scs_printf("ERROR: GPU indirect solver: out of host memory\n");
    return -1;
  }

  /* Get handle to the CUBLAS context */
  if (cublasCreate(&p->cublas_handle) != CUBLAS_STATUS_SUCCESS) {
    p->cublas_handle = 0;
    scs_printf("ERROR: GPU indirect solver: could not create the cuBLAS "
               "context\n");
    return -1;
  }
  /* Get handle to the CUSPARSE context */
  if (cusparseCreate(&p->cusparse_handle) != CUSPARSE_STATUS_SUCCESS) {
    p->cusparse_handle = 0;
    scs_printf("ERROR: GPU indirect solver: could not create the cuSPARSE "
               "context\n");
    return -1;
  }

  if (!(p->Ag = gpu_matrix(A, "A"))) {
    return -1;
  }
  if (P) {
    ScsMatrix *P_full = fill_p_matrix(P);
    if (!P_full) {
      scs_printf("ERROR: GPU indirect solver: out of host memory\n");
      return -1;
    }
    p->Pg = gpu_matrix(P_full, "P");
    SCS(cs_spfree)(P_full);
    if (!p->Pg) {
      return -1;
    }
  }

  if (gpu_malloc((void **)&p->p, n, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->r, n, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->Gp, n, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->bg, n + m, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->tmp_m, m, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->z, n, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->M_gpu, n, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->r_x_gpu, n, sizeof(scs_float), vecs) ||
      gpu_malloc((void **)&p->inv_r_y_gpu, m, sizeof(scs_float), vecs)) {
    return -1;
  }

  /* we initialize with tmp_m but always overwrite it so it doesn't matter */
  if (cusparseCreateDnVec(&p->dn_vec_n, n, p->tmp_m, SCS_CUDA_FLOAT) !=
          CUSPARSE_STATUS_SUCCESS ||
      cusparseCreateDnVec(&p->dn_vec_n_p, n, p->tmp_m, SCS_CUDA_FLOAT) !=
          CUSPARSE_STATUS_SUCCESS ||
      cusparseCreateDnVec(&p->dn_vec_m, m, p->tmp_m, SCS_CUDA_FLOAT) !=
          CUSPARSE_STATUS_SUCCESS) {
    scs_printf("ERROR: GPU indirect solver: cuSPARSE could not describe the "
               "vectors\n");
    return -1;
  }

#if GPU_TRANSPOSE_MAT > 0
  if (transpose_a(p)) {
    return -1;
  }
#endif
  if (reserve_spmv_buffer(p)) {
    return -1;
  }

  /* eigCG deflation, with all of its GPU memory (linsys/eigcg.h): W, AW and
   * the window (SCS(eigcg_init)), the diagonal for its filter and room for
   * its small matrices, (4 DEFLATE_VECTORS + window + 1) n numbers in all.
   * As on the CPU, setup fails, with a message, if that does not fit. */
  if (SCS(eigcg_init)(&p->dfl, p, &gpu_eig_ops, n) != 0) {
    return -1;
  }
  if (p->dfl.max > 0) {
    p->dfl_small_n = (size_t)p->dfl.win * 2 * p->dfl.max;
    if (gpu_malloc((void **)&p->Md_gpu, n, sizeof(scs_float),
                   "eigCG deflation's diagonal") ||
        gpu_malloc((void **)&p->dfl_small, p->dfl_small_n, sizeof(scs_float),
                   "eigCG deflation's small matrices")) {
      return -1;
    }
  }

  /* Form preconditioner and copy R_x, 1/R_y to gpu */
  if (scs_update_lin_sys_diag_r(p, diag_r) != 0) {
    return -1;
  }

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    scs_printf("%s:%d:%s\nERROR_CUDA (*): %s\n", __FILE__, __LINE__, __func__,
               cudaGetErrorString(err));
    return -1;
  }
  return SCS(lin_sys_shared_init)(&p->sh, A, P, diag_r, 1, plain_solve,
                                  SCS_NULL, SCS_NULL);
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  cudaError_t err;
  ScsLinSysWork *p = SCS_NULL;
  int device_count;

  err = cudaGetDeviceCount(&device_count);
  if (err > 0) {
    scs_printf("cudaError: %i (100 indicates no device)\n", (int)err);
    return SCS_NULL;
  }

  p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p) {
    scs_printf("ERROR: GPU indirect solver: out of host memory\n");
    return SCS_NULL;
  }
  p->A = A;
  p->P = P;
  p->m = A->m;
  p->n = A->n;
  p->tot_cg_its = 0;
  p->buffer_size = 0;
  p->buffer = SCS_NULL;

  if (init_work(p, diag_r) != 0) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  return p;
}

/* ======================== Conjugate Gradient Solver ======================== */

/* solves (R_x + P + A' R_y^{-1} A)x = b, s warm start, solution stored in
 * b, on GPU; returns the iterations, or -1 if a cuBLAS call of the deflation
 * failed */
static scs_int pcg(ScsLinSysWork *pr, const scs_float *s, scs_float *bg,
                   scs_int max_its, scs_float tol) {
  scs_int i, n = pr->n;
  scs_float ztr, ztr_prev, alpha, ptGp, beta, neg_alpha;
  scs_float onef = 1.0, neg_onef = -1.0;
  scs_float *p = pr->p;   /* cg direction */
  scs_float *Gp = pr->Gp; /* updated CG direction */
  scs_float *r = pr->r;   /* cg residual */
  scs_float *z = pr->z;   /* preconditioned */
  cublasHandle_t cublas_handle = pr->cublas_handle;

  if (!s) {
    /* take s = 0 */
    /* r = b */
    cudaMemcpy(r, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    /* b = 0 */
    cudaMemset(bg, 0, n * sizeof(scs_float));
  } else {
    /* p contains bg temporarily */
    cudaMemcpy(p, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    /* bg = s */
    cudaMemcpy(bg, s, n * sizeof(scs_float), cudaMemcpyHostToDevice);
    /* r = Mat * s */
    mat_vec(pr, bg, r);
    /* r = Mat * s - b */
    CUBLAS(axpy)(cublas_handle, n, &neg_onef, p, 1, r, 1);
    /* r = b - Mat * s  */
    CUBLAS(scal)(cublas_handle, n, &neg_onef, r, 1);
  }

  /* deflation: project the residual against the harvested eigenvectors
   * (linsys/eigcg.h) */
  if (SCS(eigcg_deflate)(&pr->dfl, bg, r) != 0) {
    return -1;
  }

  /* check to see if we need to run CG at all */
  if (cg_gpu_norm(cublas_handle, r, n) < tol) {
    return 0;
  }

  /* z = M r */
  cudaMemcpy(z, r, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
  scale_by_diag(cublas_handle, pr->M_gpu, z, n);
  /* ztr = z'r */
  CUBLAS(dot)(cublas_handle, n, r, 1, z, 1, &ztr);
  /* p = z */
  cudaMemcpy(p, z, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
  SCS(eigcg_start)(&pr->dfl, z, ztr);

  for (i = 0; i < max_its; ++i) {
    /* Gp = Mat * p */
    mat_vec(pr, p, Gp);
    /* ptGp = p'Gp */
    CUBLAS(dot)(cublas_handle, n, p, 1, Gp, 1, &ptGp);
    /* alpha = z'r / p'G p */
    alpha = ztr / ptGp;
    neg_alpha = -alpha;
    /* b += alpha * p */
    CUBLAS(axpy)(cublas_handle, n, &alpha, p, 1, bg, 1);
    /* r -= alpha * G p */
    CUBLAS(axpy)(cublas_handle, n, &neg_alpha, Gp, 1, r, 1);

#if VERBOSITY > 3
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol,
               cg_gpu_norm(cublas_handle, r, n), (long)i + 1);
#endif

    if (cg_gpu_norm(cublas_handle, r, n) < tol) {
      return i + 1;
    }
    /* z = M r */
    cudaMemcpy(z, r, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    scale_by_diag(cublas_handle, pr->M_gpu, z, n);
    ztr_prev = ztr;
    /* ztr = z'r */
    CUBLAS(dot)(cublas_handle, n, r, 1, z, 1, &ztr);
    if (ztr_prev == 0.) {
      break;
    }
    beta = ztr / ztr_prev;
    SCS(eigcg_step)(&pr->dfl, alpha, beta, z, ztr);
    /* p = beta * p, where beta = ztr / ztr_prev */
    CUBLAS(scal)(cublas_handle, n, &beta, p, 1);
    /* p = z + beta * p */
    CUBLAS(axpy)(cublas_handle, n, &onef, z, 1, p, 1);
  }
  return i;
}

/* solves Mx = b, for x but stores result in b */
/* s contains warm-start (if available) */
/*
 * [x] = [R_x + P        A' ]^{-1} [rx]
 * [y]   [     A        -R_y ]      [ry]
 *
 * becomes:
 *
 * x = (R_x + P + A' R_y^{-1} A)^{-1} (rx + A' R_y^{-1} ry)
 * y = R_y^{-1} (Ax - ry)
 *
 */
static scs_int plain_solve(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol) {
  scs_int cg_its, max_iters, probes;
  scs_float neg_onef = -1.0;

  /* these are on GPU */
  scs_float *bg = p->bg;
  scs_float *tmp_m = p->tmp_m;
  ScsGpuMatrix *Ag = p->Ag;

  /* the cold solves (s == SCS_NULL: the solve for g, and the refinement
   * steps of scs_solve_lin_sys) are the ones worth harvesting from, the
   * first one only (no set held); warm-started solves do not harvest */
  SCS(eigcg_begin_solve)(&p->dfl, !s);

  /* x = 0 if b is this small, see CG_ZERO_RHS_TOL */
  if (CG_NORM(b, p->n + p->m) <= MIN(tol, CG_ZERO_RHS_TOL)) {
    memset(b, 0, (p->n + p->m) * sizeof(scs_float));
    return 0;
  }

  if (tol <= 0.) {
    scs_printf("Warning: tol = %4f <= 0, likely compiled without setting "
               "INDIRECT flag.\n",
               tol);
  }

  /* bg = b = [rx; ry] */
  cudaMemcpy(bg, b, (Ag->n + Ag->m) * sizeof(scs_float),
             cudaMemcpyHostToDevice);
  /* tmp = ry */
  cudaMemcpy(tmp_m, &(bg[Ag->n]), Ag->m * sizeof(scs_float),
             cudaMemcpyDeviceToDevice);
  /* tmp = R_y^{-1} * tmp = R_y^{-1} * ry */
  scale_by_diag(p->cublas_handle, p->inv_r_y_gpu, tmp_m, p->Ag->m);

  cusparseDnVecSetValues(p->dn_vec_m, (void *)tmp_m); /* R * ry */
  cusparseDnVecSetValues(p->dn_vec_n, (void *)bg);    /* rx */
  /* bg[:n] = rx + A' R ry */
  SCS(accum_by_atrans_gpu)(Ag, p->dn_vec_m, p->dn_vec_n, p->cusparse_handle,
                           &p->buffer_size, &p->buffer);

  /* set max_iters to 10 * n (though in theory n is enough for any tol) */
  max_iters = 10 * Ag->n;

  /* solves (R_x + P + A' R_y^{-1} A)x = bg, s warm start, solution stored
   * in bg */
  cg_its = pcg(p, s, bg, max_iters, tol); /* bg[:n] = x */
  /* if it harvested, extract the Ritz vectors; the image matvecs it spends
   * are added to the iteration count so the matvec accounting stays honest.
   * It uses CG's vectors and the descriptors as scratch, and leaves bg. */
  probes = SCS(eigcg_end_solve)(&p->dfl);
  if (cg_its < 0 || probes < 0) {
    scs_printf("ERROR: GPU indirect solver: a cuBLAS call of eigCG "
               "deflation failed\n");
    return -1;
  }
  cg_its += probes;

  /* bg[n:] = -ry */
  CUBLAS(scal)(p->cublas_handle, Ag->m, &neg_onef, &(bg[Ag->n]), 1);
  cusparseDnVecSetValues(p->dn_vec_m, (void *)&(bg[Ag->n])); /* -ry */
  cusparseDnVecSetValues(p->dn_vec_n, (void *)bg);           /* x */

  /* b[n:] = Ax - ry */
#if GPU_TRANSPOSE_MAT > 0
  SCS(accum_by_atrans_gpu)(p->Agt, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle,
                           &p->buffer_size, &p->buffer);
#else
  SCS(accum_by_a_gpu)(Ag, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle,
                      &p->buffer_size, &p->buffer);
#endif

  /* bg[n:] = R_y^{-1} bg[n:] = R_y^{-1} (Ax - ry) = y */
  scale_by_diag(p->cublas_handle, p->inv_r_y_gpu, &(bg[p->n]), p->Ag->m);

  /* copy bg = [x; y] back to b */
  cudaMemcpy(b, bg, (Ag->n + Ag->m) * sizeof(scs_float),
             cudaMemcpyDeviceToHost);
  p->tot_cg_its += cg_its;
#if VERBOSITY > 1
  scs_printf("tol %.3e\n", tol);
  scs_printf("cg_its %i\n", (int)cg_its);
#endif
  return 0;
}

/* the warm start, refinement and residual around plain_solve, see linsys.h;
 * the residuals are SCS's products on the CPU */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol, scs_int max_refine_steps,
                          scs_float *res) {
  return SCS(lin_sys_shared_solve)(&p->sh, p, b, s, tol, max_refine_steps,
                                   res);
}
