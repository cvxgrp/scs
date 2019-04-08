#include "private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

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
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str, "sparse-indirect GPU, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)",
          (long)A->p[A->n], stgs->cg_rate);
  return str;
}

char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str,
          "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
          (scs_float)p->tot_cg_its / (info->iter + 1),
          p->total_solve_time / (info->iter + 1) / 1e3);
  p->tot_cg_its = 0;
  p->total_solve_time = 0;
  return str;
}

static void cuda_free_a_matrix(ScsMatrix *A) {
  if (A->x) {
    cudaFree(A->x);
  }
  if (A->i) {
    cudaFree(A->i);
  }
  if (A->p) {
    cudaFree(A->p);
  }
}

void SCS(free_lin_sys_work)(ScsLinSysWork *p) {
  if (p) {
    if (p->p) {
      cudaFree(p->p);
    }
    if (p->r) {
      cudaFree(p->r);
    }
    if (p->Gp) {
      cudaFree(p->Gp);
    }
    if (p->bg) {
      cudaFree(p->bg);
    }
    if (p->tmp_m) {
      cudaFree(p->tmp_m);
    }
    if (p->z) {
      cudaFree(p->z);
    }
    if (p->M) {
      cudaFree(p->M);
    }
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
    /* Don't reset because it interferes with other GPU programs. */
    /* cudaDeviceReset(); */
    scs_free(p);
  }
}

/*y = (RHO_X * I + A'A)x */
static void mat_vec(const ScsMatrix *A, const ScsSettings *s, ScsLinSysWork *p,
                    const scs_float *x, scs_float *y) {
  /* x and y MUST already be loaded to GPU */
  scs_float *tmp_m = p->tmp_m; /* temp memory */
  cudaMemset(tmp_m, 0, A->m * sizeof(scs_float));
  accum_by_a_gpu(p, x, tmp_m);
  cudaMemset(y, 0, A->n * sizeof(scs_float));
  accum_by_atrans_gpu(p, tmp_m, y);
  CUBLAS(axpy)(p->cublas_handle, A->n, &(s->rho_x), x, 1, y, 1);
}

/* M = inv ( diag ( RHO_X * I + A'A ) ) */
static void get_preconditioner(const ScsMatrix *A, const ScsSettings *stgs,
                               ScsLinSysWork *p) {
  scs_int i;
  scs_float *M = (scs_float *)scs_malloc(A->n * sizeof(scs_float));

#if EXTRA_VERBOSE > 0
  scs_printf("getting pre-conditioner\n");
#endif

  for (i = 0; i < A->n; ++i) {
    M[i] = 1 / (stgs->rho_x +
                SCS(norm_sq)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
    /* M[i] = 1; */
  }
  cudaMemcpy(p->M, M, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
  scs_free(M);

#if EXTRA_VERBOSE > 0
  scs_printf("finished getting pre-conditioner\n");
#endif
}

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
  p->tot_cg_its = 0;

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

  cudaMalloc((void **)&p->p, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->r, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->Gp, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->bg, (A->n + A->m) * sizeof(scs_float));
  cudaMalloc((void **)&p->tmp_m,
             A->m * sizeof(scs_float)); /* intermediate result */
  cudaMalloc((void **)&p->z, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->M, A->n * sizeof(scs_float));

  cudaMemcpy(Ag->i, A->i, (A->p[A->n]) * sizeof(scs_int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(Ag->p, A->p, (A->n + 1) * sizeof(scs_int), cudaMemcpyHostToDevice);
  cudaMemcpy(Ag->x, A->x, (A->p[A->n]) * sizeof(scs_float),
             cudaMemcpyHostToDevice);

  get_preconditioner(A, stgs, p);

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

static void apply_pre_conditioner(cublasHandle_t cublas_handle, scs_float *M,
                                  scs_float *z, scs_float *r, scs_int n) {
  cudaMemcpy(z, r, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
  CUBLAS(tbmv)
  (cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n,
   0, M, 1, z, 1);
}

/* solves (I+A'A)x = b, s warm start, solution stored in bg (on GPU) */
static scs_int pcg(const ScsMatrix *A, const ScsSettings *stgs,
                   ScsLinSysWork *pr, const scs_float *s, scs_float *bg,
                   scs_int max_its, scs_float tol) {
  scs_int i, n = A->n;
  scs_float alpha, nrm_r, p_gp, neg_alpha, beta, ipzr, ipzr_old;
  scs_float onef = 1.0, neg_onef = -1.0;
  scs_float *p = pr->p;   /* cg direction */
  scs_float *Gp = pr->Gp; /* updated CG direction */
  scs_float *r = pr->r;   /* cg residual */
  scs_float *z = pr->z;   /* preconditioned */
  scs_float *M = pr->M;   /* preconditioner */
  cublasHandle_t cublas_handle = pr->cublas_handle;

  if (s == SCS_NULL) {
    cudaMemcpy(r, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    cudaMemset(bg, 0, n * sizeof(scs_float));
  } else {
    /* p contains bg temporarily */
    cudaMemcpy(p, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    /* bg contains s */
    cudaMemcpy(bg, s, n * sizeof(scs_float), cudaMemcpyHostToDevice);
    mat_vec(A, stgs, pr, bg, r);
    CUBLAS(axpy)(cublas_handle, n, &neg_onef, p, 1, r, 1);
    CUBLAS(scal)(cublas_handle, n, &neg_onef, r, 1);
  }

  /* for some reason nrm2 is VERY slow */
  /* CUBLAS(nrm2)(cublas_handle, n, r, 1, &nrm_r); */
  CUBLAS(dot)(cublas_handle, n, r, 1, r, 1, &nrm_r);
  nrm_r = SQRTF(nrm_r);
  /* check to see if we need to run CG at all */
  if (nrm_r < MIN(tol, 1e-18)) {
    return 0;
  }

  apply_pre_conditioner(cublas_handle, M, z, r, n);
  CUBLAS(dot)(cublas_handle, n, r, 1, z, 1, &ipzr);
  /* put z in p, replacing temp mem */
  cudaMemcpy(p, z, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);

  for (i = 0; i < max_its; ++i) {
    mat_vec(A, stgs, pr, p, Gp);

    CUBLAS(dot)(cublas_handle, n, p, 1, Gp, 1, &p_gp);

    alpha = ipzr / p_gp;
    neg_alpha = -alpha;

    CUBLAS(axpy)(cublas_handle, n, &alpha, p, 1, bg, 1);
    CUBLAS(axpy)(cublas_handle, n, &neg_alpha, Gp, 1, r, 1);

    /* for some reason nrm2 is VERY slow */
    /* CUBLAS(nrm2)(cublas_handle, n, r, 1, &nrm_r); */
    CUBLAS(dot)(cublas_handle, n, r, 1, r, 1, &nrm_r);
    nrm_r = SQRTF(nrm_r);
    if (nrm_r < tol) {
      i++;
      break;
    }
    ipzr_old = ipzr;
    apply_pre_conditioner(cublas_handle, M, z, r, n);
    CUBLAS(dot)(cublas_handle, n, r, 1, z, 1, &ipzr);

    beta = ipzr / ipzr_old;
    CUBLAS(scal)(cublas_handle, n, &beta, p, 1);
    CUBLAS(axpy)(cublas_handle, n, &onef, z, 1, p, 1);
  }
#if EXTRA_VERBOSE > 0
  scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, nrm_r, (long)i + 1);
#endif
  return i;
}

#ifdef TEST_GPU_MAT_MUL
static void accum_by_atrans_host(const ScsMatrix *A, ScsLinSysWork *p,
                                 const scs_float *x, scs_float *y) {
  SCS(_accum_by_atrans)(A->n, A->x, A->i, A->p, x, y);
}

static void accum_by_a_host(const ScsMatrix *A, ScsLinSysWork *p,
                            const scs_float *x, scs_float *y) {
  SCS(_accum_by_a)(A->n, A->x, A->i, A->p, x, y);
}

static void test_gpu_mat_mul(const ScsMatrix *A, ScsLinSysWork *p,
                             scs_float *b) {
  /* test to see if matrix multiplication codes agree */
  scs_float t[A->n + A->m], u[A->n + A->m], *bg;
  cudaMalloc((void **)&bg, (A->n + A->m) * sizeof(scs_float));

  cudaMemcpy(bg, b, (A->n + A->m) * sizeof(scs_float), cudaMemcpyHostToDevice);
  memcpy(t, b, (A->n + A->m) * sizeof(scs_float));

  accum_by_atrans_gpu(p, &(bg[A->n]), bg);
  accum_by_atrans_host(A, p, &(t[A->n]), t);
  cudaMemcpy(u, bg, (A->n + A->m) * sizeof(scs_float), cudaMemcpyDeviceToHost);
  scs_printf("A trans multiplication err %2.e\n", SCS(norm_diff)(u, t, A->n));

  accum_by_a_gpu(p, bg, &(bg[A->n]));
  accum_by_a_host(A, p, t, &(t[A->n]));
  cudaMemcpy(u, bg, (A->n + A->m) * sizeof(scs_float), cudaMemcpyDeviceToHost);
  scs_printf("A multiplcation err %2.e\n",
             SCS(norm_diff)(&(u[A->n]), &(t[A->n]), A->m));
  cudaFree(bg);
}
#endif

scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsSettings *stgs,
                           ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_int iter) {
  scs_int cg_its;
  SCS(timer) linsys_timer;
  scs_float *bg = p->bg;
  scs_float neg_onef = -1.0;
  scs_float cg_tol =
      SCS(norm)(b, A->n) *
      (iter < 0 ? CG_BEST_TOL
                : CG_MIN_TOL / POWF((scs_float)iter + 1., stgs->cg_rate));

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
