#include "private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CUDA_CHECK_ERR                                                    \
  do {                                                                    \
    cudaError_t err = cudaGetLastError();                                 \
    if (err != cudaSuccess) {                                             \
      printf("%s:%d:%s\n ERROR_CUDA: %s\n", __FILE__, __LINE__, __func__, \
             cudaGetErrorString(err));                                    \
    }                                                                     \
  } while (0)

#ifndef EXTRA_VERBOSE
#ifndef SFLOAT
#define CUBLAS(x) cublasD##x
#define CUSPARSE(x) cusparseD##x
#else
#define CUBLAS(x) cublasS##x
#define CUSPARSE(x) cusparseS##x
#endif
#else
#ifndef SFLOAT
#define CUBLAS(x) \
  CUDA_CHECK_ERR; \
  cublasD##x
#define CUSPARSE(x) \
  CUDA_CHECK_ERR;   \
  cusparseD##x
#else
#define CUBLAS(x) \
  CUDA_CHECK_ERR; \
  cublasS##x
#define CUSPARSE(x) \
  CUDA_CHECK_ERR;   \
  cusparseS##x
#endif
#endif

/*
 CUDA matrix routines only for CSR, not CSC matrices:
    CSC             CSR             GPU     Mult
    A (m x n)       A' (n x m)      Ag      accum_by_a_trans_gpu
    A'(n x m)       A  (m x n)      Agt     accum_by_a_gpu
*/

static void accum_by_atrans_gpu(const ScsLinSysWork *p, const scs_float *x,
                                scs_float *y) {
  /* y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  ScsMatrix *Ag = p->Ag;
  CUSPARSE(csrmv)
  (p->cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Ag->n, Ag->m, p->Annz,
   &onef, p->descr, Ag->x, Ag->p, Ag->i, x, &onef, y);
}

static void accum_by_a_gpu(const ScsLinSysWork *p, const scs_float *x,
                           scs_float *y) {
  /* y += A*x
     x and y MUST be on GPU already
   */
  const scs_float onef = 1.0;
#if GPU_TRANSPOSE_MAT > 0
  ScsMatrix *Agt = p->Agt;
  CUSPARSE(csrmv)
  (p->cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Agt->n, Agt->m,
   p->Annz, &onef, p->descr, Agt->x, Agt->p, Agt->i, x, &onef, y);
#else
  /* The A matrix idx pointers must be ORDERED */
  ScsMatrix *Ag = p->Ag;
  CUSPARSE(csrmv)
  (p->cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE, Ag->n, Ag->m, p->Annz,
   &onef, p->descr, Ag->x, Ag->p, Ag->i, x, &onef, y);
#endif
}

/* do not use within pcg, reuses memory */
void SCS(accum_by_atrans)(const ScsMatrix *A, ScsLinSysWork *p,
                          const scs_float *x, scs_float *y) {
  scs_float *v_m = p->tmp_m;
  scs_float *v_n = p->r;
  cudaMemcpy(v_m, x, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
  cudaMemcpy(v_n, y, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
  accum_by_atrans_gpu(p, v_m, v_n);
  cudaMemcpy(y, v_n, A->n * sizeof(scs_float), cudaMemcpyDeviceToHost);
}

/* do not use within pcg, reuses memory */
void SCS(accum_by_a)(const ScsMatrix *A, ScsLinSysWork *p, const scs_float *x,
                     scs_float *y) {
  scs_float *v_m = p->tmp_m;
  scs_float *v_n = p->r;
  cudaMemcpy(v_n, x, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
  cudaMemcpy(v_m, y, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
  accum_by_a_gpu(p, v_n, v_m);
  cudaMemcpy(y, v_m, A->m * sizeof(scs_float), cudaMemcpyDeviceToHost);
}

char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsSettings *stgs) {
  char *tmp = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(tmp, "sparse-direct GPU, nnz in A = %li", (long)A->p[A->n]);
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

static void cuda_free_a_matrix(ScsMatrix *A) {
    cudaFree(A->x);
    cudaFree(A->i);
    cudaFree(A->p);
}

// TODO
void SCS(free_lin_sys_work)(ScsLinSysWork *p) {
  if (p) {
      cudaFree(p->p);
      cudaFree(p->r);
      cudaFree(p->Gp);
      cudaFree(p->bg);
      cudaFree(p->tmp_m);
      cudaFree(p->z);
      cudaFree(p->M);
    if (p->Ag) {
      cuda_free_a_matrix(p->Ag);
      scs_free(p->Ag);
    }
    if (p->Agt) {
      cuda_free_a_matrix(p->Agt);
      scs_free(p->Agt);
    }
    cusparseDestroy(p->cusparse_handle);
    cublasDestroy(p->cublas_handle);
    cusparseDestroyMatDescr(p->descr);
    /* Don't reset because it interferes with other GPU programs. */
    /* cudaDeviceReset(); */
    scs_free(p);
  }
}
/* This is copied directly from SCS direct CPU code. */
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

/* This is copied directly from SCS direct CPU code. */
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

/* This is copied directly from SCS direct CPU code. */
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

// TODO
ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A,
                                      const ScsSettings *stgs) {
  cudaError_t err;
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  ScsMatrix *Ag = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));

  p->Annz = A->p[A->n];
  p->cublas_handle = 0;
  p->cusparse_handle = 0;
  p->descr = 0;

  p->total_solve_time = 0;

  /* Get handle to the CUBLAS context */
  cublasCreate(&p->cublas_handle);

  /* Get handle to the CUSPARSE context */
  cusparseCreate(&p->cusparse_handle);

  /* Matrix description */
  cusparseCreateMatDescr(&p->descr);
  cusparseSetMatType(p->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(p->descr, CUSPARSE_INDEX_BASE_ZERO);

  Ag->n = A->n;
  Ag->m = A->m;
  p->Ag = Ag;
  p->Agt = SCS_NULL;

  cudaMalloc((void **)&Ag->i, (A->p[A->n]) * sizeof(scs_int));
  cudaMalloc((void **)&Ag->p, (A->n + 1) * sizeof(scs_int));
  cudaMalloc((void **)&Ag->x, (A->p[A->n]) * sizeof(scs_float));

  cudaMemcpy(Ag->i, A->i, (A->p[A->n]) * sizeof(scs_int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(Ag->p, A->p, (A->n + 1) * sizeof(scs_int), cudaMemcpyHostToDevice);
  cudaMemcpy(Ag->x, A->x, (A->p[A->n]) * sizeof(scs_float),
             cudaMemcpyHostToDevice);

#if GPU_TRANSPOSE_MAT > 0
  p->Agt = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));
  p->Agt->n = A->m;
  p->Agt->m = A->n;
  cudaMalloc((void **)&p->Agt->i, (A->p[A->n]) * sizeof(scs_int));
  cudaMalloc((void **)&p->Agt->p, (A->m + 1) * sizeof(scs_int));
  cudaMalloc((void **)&p->Agt->x, (A->p[A->n]) * sizeof(scs_float));
  /* transpose Ag into Agt for faster multiplies */
  /* TODO: memory intensive, could perform transpose in CPU and copy to GPU */
  CUSPARSE(csr2csc)
  (p->cusparse_handle, A->n, A->m, A->p[A->n], Ag->x, Ag->p, Ag->i, p->Agt->x,
   p->Agt->i, p->Agt->p, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
#endif

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("%s:%d:%s\nERROR_CUDA: %s\n", __FILE__, __LINE__, __func__,
           cudaGetErrorString(err));
    SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  return p;
}

scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsSettings *stgs,
                           ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_int iter) {
  scs_int cg_its;
  SCS(timer) linsys_timer;
  scs_float *bg = p->bg;
  scs_float neg_onef = -1.0;
  SCS(tic)(&linsys_timer);

  /* solves Mx = b, for x but stores result in b */
  /* s contains warm-start (if available) */

#ifdef TEST_GPU_MAT_MUL
  test_gpu_mat_mul(A, p, b);
#endif

  /* all on GPU */
  cudaMemcpy(bg, b, (A->n + A->m) * sizeof(scs_float), cudaMemcpyHostToDevice);
  accum_by_atrans_gpu(p, &(bg[A->n]), bg);
  /* solves (I+A'A)x = b, s warm start, solution stored in b */
  cg_its = pcg(A, stgs, p, s, bg, A->n, MAX(cg_tol, CG_BEST_TOL));
  CUBLAS(scal)(p->cublas_handle, A->m, &neg_onef, &(bg[A->n]), 1);
  accum_by_a_gpu(p, bg, &(bg[A->n]));
  cudaMemcpy(b, bg, (A->n + A->m) * sizeof(scs_float), cudaMemcpyDeviceToHost);

  if (iter >= 0) {
    p->tot_cg_its += cg_its;
  }

  p->total_solve_time += SCS(tocq)(&linsys_timer);
#if EXTRAVERBOSE > 0
  scs_printf("linsys solve time: %1.2es\n", SCS(tocq)(&linsys_timer) / 1e3);
#endif
  return 0;
}

void SCS(normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                      const ScsCone *k, ScsScaling *scal) {
  SCS(_normalize_a)(A, stgs, k, scal);
}

void SCS(un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                         const ScsScaling *scal) {
  SCS(_un_normalize_a)(A, stgs, scal);
}

#ifdef __cplusplus
}
#endif
