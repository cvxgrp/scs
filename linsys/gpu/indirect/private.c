#include "private.h"
#include "linsys.h"

/* norm to use when deciding convergence */
/* should be consistent with CG_NORM in glbopts.h */
#define USE_L2_NORM (0)

static scs_float cg_gpu_norm(cublasHandle_t cublas_handle, scs_float *r,
                             scs_int n) {
#if USE_L2_NORM > 0
  scs_float nrm;
  CUBLAS(nrm2)(cublas_handle, n, r, 1, &nrm);
#else
  scs_int idx;
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
/* set M = inv ( diag ( R_x + P + A' R_y^{-1} A ) ) */
static void set_preconditioner(ScsLinSysWork *p, const scs_float *diag_r) {
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
    /* finally invert for pre-conditioner */
    M[i] = 1. / M[i];
  }
  cudaMemcpy(p->M_gpu, M, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
#if VERBOSITY > 0
  scs_printf("finished getting pre-conditioner\n");
#endif
}

/* no need to update anything in this case */
void scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;

  /* R_x to gpu */
  cudaMemcpy(p->r_x_gpu, diag_r, p->n * sizeof(scs_float),
             cudaMemcpyHostToDevice);

  /* 1/R_y to gpu */
  for (i = 0; i < p->m; ++i)
    p->inv_r_y[i] = 1. / diag_r[p->n + i];
  cudaMemcpy(p->inv_r_y_gpu, p->inv_r_y, p->m * sizeof(scs_float),
             cudaMemcpyHostToDevice);

  /* set preconditioner M on gpu */
  set_preconditioner(p, diag_r);
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    scs_free(p->M);
    scs_free(p->inv_r_y);
    cudaFree(p->p);
    cudaFree(p->r);
    cudaFree(p->Gp);
    cudaFree(p->bg);
    cudaFree(p->tmp_m);
    cudaFree(p->z);
    cudaFree(p->M_gpu);
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
    cusparseDestroyDnVec(p->dn_vec_m);
    cusparseDestroyDnVec(p->dn_vec_n);
    cusparseDestroyDnVec(p->dn_vec_n_p);
    cusparseDestroy(p->cusparse_handle);
    cublasDestroy(p->cublas_handle);
    /* Don't reset because it interferes with other GPU programs. */
    /* cudaDeviceReset(); */
    scs_free(p);
  }
}

/* z = M * z elementwise in place, assumes M, z on GPU */
static void scale_by_diag(cublasHandle_t cublas_handle, scs_float *M,
                          scs_float *z, scs_int n) {
  CUBLAS(tbmv)
  (cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n,
   0, M, 1, z, 1);
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
  cudaMemcpy(y, x, p->n * sizeof(scs_float), cudaMemcpyHostToDevice);
  /* y = R_x * x */
  scale_by_diag(p->cublas_handle, p->r_x_gpu, y, p->n);

  if (p->Pg) {
    /* y = R_x * x + P x */
    SCS(accum_by_p_gpu)
    (p->Pg, p->dn_vec_n, p->dn_vec_n_p, p->cusparse_handle, &p->buffer_size,
     &p->buffer);
  }

  /* z = Ax */
#if GPU_TRANSPOSE_MAT > 0
  SCS(accum_by_atrans_gpu)
  (p->Agt, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle, &p->buffer_size,
   &p->buffer);
#else
  SCS(accum_by_a_gpu)
  (p->Ag, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle, &p->buffer_size,
   &p->buffer);
#endif
  /* z = R_y^{-1} A x */
  scale_by_diag(p->cublas_handle, p->inv_r_y_gpu, z, p->m);

  /* y += A'z => y = R_x * x + P x + A' R_y^{-1} Ax */
  SCS(accum_by_atrans_gpu)
  (p->Ag, p->dn_vec_m, p->dn_vec_n_p, p->cusparse_handle, &p->buffer_size,
   &p->buffer);
}

/* P comes in upper triangular, expand to full
 * First compute triplet version of full matrix, then compress to CSC
 * */
static ScsMatrix *fill_p_matrix(const ScsMatrix *P) {
  scs_int i, j, k, kk;
  scs_int Pnzmax = 2 * P->p[P->n]; /* upper bound */
  ScsMatrix *P_tmp = SCS(cs_spalloc)(P->n, P->n, Pnzmax, 1, 1);
  ScsMatrix *P_full;
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

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  cudaError_t err;
  ScsMatrix *P_full;
  ScsLinSysWork *p = SCS_NULL;
  ScsGpuMatrix *Ag = SCS_NULL;
  ScsGpuMatrix *Pg = SCS_NULL;
  int device_count;

  err = cudaGetDeviceCount(&device_count);
  if (err > 0) {
    scs_printf("cudaError: %i (100 indicates no device)\n", (int)err);
    return SCS_NULL;
  }

  p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  Ag = (ScsGpuMatrix *)scs_calloc(1, sizeof(ScsGpuMatrix));

  p->inv_r_y = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  p->M = (scs_float *)scs_calloc(A->n, sizeof(scs_float));

  p->A = A;
  p->P = P;
  p->m = A->m;
  p->n = A->n;

#if GPU_TRANSPOSE_MAT > 0
  size_t new_buffer_size = 0;
#endif

  p->cublas_handle = 0;
  p->cusparse_handle = 0;

  p->tot_cg_its = 0;

  p->buffer_size = 0;
  p->buffer = SCS_NULL;

  /* Get handle to the CUBLAS context */
  cublasCreate(&p->cublas_handle);

  /* Get handle to the CUSPARSE context */
  cusparseCreate(&p->cusparse_handle);

  Ag->n = A->n;
  Ag->m = A->m;
  Ag->nnz = A->p[A->n];
  Ag->descr = 0;
  cudaMalloc((void **)&Ag->i, (A->p[A->n]) * sizeof(scs_int));
  cudaMalloc((void **)&Ag->p, (A->n + 1) * sizeof(scs_int));
  cudaMalloc((void **)&Ag->x, (A->p[A->n]) * sizeof(scs_float));

  cudaMalloc((void **)&p->p, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->r, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->Gp, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->bg, (A->n + A->m) * sizeof(scs_float));
  cudaMalloc((void **)&p->tmp_m, A->m * sizeof(scs_float));
  cudaMalloc((void **)&p->z, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->M_gpu, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->r_x_gpu, A->n * sizeof(scs_float));
  cudaMalloc((void **)&p->inv_r_y_gpu, A->m * sizeof(scs_float));

  cudaMemcpy(Ag->i, A->i, (A->p[A->n]) * sizeof(scs_int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(Ag->p, A->p, (A->n + 1) * sizeof(scs_int), cudaMemcpyHostToDevice);
  cudaMemcpy(Ag->x, A->x, (A->p[A->n]) * sizeof(scs_float),
             cudaMemcpyHostToDevice);

  cusparseCreateCsr(&Ag->descr, Ag->n, Ag->m, Ag->nnz, Ag->p, Ag->i, Ag->x,
                    SCS_CUSPARSE_INDEX, SCS_CUSPARSE_INDEX,
                    CUSPARSE_INDEX_BASE_ZERO, SCS_CUDA_FLOAT);

  if (P) {
    Pg = (ScsGpuMatrix *)scs_calloc(1, sizeof(ScsGpuMatrix));
    P_full = fill_p_matrix(P);
    Pg->n = P_full->n;
    Pg->m = P_full->m;
    Pg->nnz = P_full->p[P_full->n];
    Pg->descr = 0;
    cudaMalloc((void **)&Pg->i, (P_full->p[P_full->n]) * sizeof(scs_int));
    cudaMalloc((void **)&Pg->p, (P_full->n + 1) * sizeof(scs_int));
    cudaMalloc((void **)&Pg->x, (P_full->p[P_full->n]) * sizeof(scs_float));

    cudaMemcpy(Pg->i, P_full->i, (P_full->p[P_full->n]) * sizeof(scs_int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(Pg->p, P_full->p, (P_full->n + 1) * sizeof(scs_int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(Pg->x, P_full->x, (P_full->p[P_full->n]) * sizeof(scs_float),
               cudaMemcpyHostToDevice);

    cusparseCreateCsr(&Pg->descr, Pg->n, Pg->m, Pg->nnz, Pg->p, Pg->i, Pg->x,
                      SCS_CUSPARSE_INDEX, SCS_CUSPARSE_INDEX,
                      CUSPARSE_INDEX_BASE_ZERO, SCS_CUDA_FLOAT);

    SCS(cs_spfree)(P_full);
  } else {
    Pg = SCS_NULL;
  }

  p->Ag = Ag;
  p->Pg = Pg;
  p->Agt = SCS_NULL;

  /* we initialize with tmp_m but always overwrite it so it doesn't matter */
  cusparseCreateDnVec(&p->dn_vec_n, Ag->n, p->tmp_m, SCS_CUDA_FLOAT);
  cusparseCreateDnVec(&p->dn_vec_n_p, Ag->n, p->tmp_m, SCS_CUDA_FLOAT);
  cusparseCreateDnVec(&p->dn_vec_m, Ag->m, p->tmp_m, SCS_CUDA_FLOAT);

  /* Form preconditioner and copy R_x, 1/R_y to gpu */
  scs_update_lin_sys_diag_r(p, diag_r);

#if GPU_TRANSPOSE_MAT > 0
  p->Agt = (ScsGpuMatrix *)scs_malloc(sizeof(ScsGpuMatrix));
  p->Agt->n = A->m;
  p->Agt->m = A->n;
  p->Agt->nnz = A->p[A->n];
  p->Agt->descr = 0;
  /* Matrix description */

  cudaMalloc((void **)&p->Agt->i, (A->p[A->n]) * sizeof(scs_int));
  cudaMalloc((void **)&p->Agt->p, (A->m + 1) * sizeof(scs_int));
  cudaMalloc((void **)&p->Agt->x, (A->p[A->n]) * sizeof(scs_float));
  /* transpose Ag into Agt for faster multiplies */
  /* TODO: memory intensive, could perform transpose in CPU and copy to GPU */
  cusparseCsr2cscEx2_bufferSize(
      p->cusparse_handle, A->n, A->m, A->p[A->n], Ag->x, Ag->p, Ag->i,
      p->Agt->x, p->Agt->p, p->Agt->i, SCS_CUDA_FLOAT, CUSPARSE_ACTION_NUMERIC,
      CUSPARSE_INDEX_BASE_ZERO, SCS_CSR2CSC_ALG, &new_buffer_size);

  if (new_buffer_size > p->buffer_size) {
    if (p->buffer != SCS_NULL) {
      cudaFree(p->buffer);
    }
    cudaMalloc(&p->buffer, new_buffer_size);
    p->buffer_size = new_buffer_size;
  }

  cusparseCsr2cscEx2(p->cusparse_handle, A->n, A->m, A->p[A->n], Ag->x, Ag->p,
                     Ag->i, p->Agt->x, p->Agt->p, p->Agt->i, SCS_CUDA_FLOAT,
                     CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO,
                     SCS_CSR2CSC_ALG, p->buffer);

  cusparseCreateCsr(&p->Agt->descr, p->Agt->n, p->Agt->m, p->Agt->nnz,
                    p->Agt->p, p->Agt->i, p->Agt->x, SCS_CUSPARSE_INDEX,
                    SCS_CUSPARSE_INDEX, CUSPARSE_INDEX_BASE_ZERO,
                    SCS_CUDA_FLOAT);
#endif

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("%s:%d:%s\nERROR_CUDA (*): %s\n", __FILE__, __LINE__, __func__,
           cudaGetErrorString(err));
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  return p;
}

/* solves (R_x + P + A' R_y^{-1} A)x = b, s warm start, solution stored in
 * b, on GPU */
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
    beta = ztr / ztr_prev;
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
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  scs_int cg_its, max_iters;
  scs_float neg_onef = -1.0;

  /* these are on GPU */
  scs_float *bg = p->bg;
  scs_float *tmp_m = p->tmp_m;
  ScsGpuMatrix *Ag = p->Ag;

  if (CG_NORM(b, p->n + p->m) <= 1e-12) {
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
  SCS(accum_by_atrans_gpu)
  (Ag, p->dn_vec_m, p->dn_vec_n, p->cusparse_handle, &p->buffer_size,
   &p->buffer);

  /* set max_iters to 10 * n (though in theory n is enough for any tol) */
  max_iters = 10 * Ag->n;

  /* solves (R_x + P + A' R_y^{-1} A)x = bg, s warm start, solution stored
   * in bg */
  cg_its = pcg(p, s, bg, max_iters, tol); /* bg[:n] = x */

  /* bg[n:] = -ry */
  CUBLAS(scal)(p->cublas_handle, Ag->m, &neg_onef, &(bg[Ag->n]), 1);
  cusparseDnVecSetValues(p->dn_vec_m, (void *)&(bg[Ag->n])); /* -ry */
  cusparseDnVecSetValues(p->dn_vec_n, (void *)bg);           /* x */

  /* b[n:] = Ax - ry */
#if GPU_TRANSPOSE_MAT > 0
  SCS(accum_by_atrans_gpu)
  (p->Agt, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle, &p->buffer_size,
   &p->buffer);
#else
  SCS(accum_by_a_gpu)
  (Ag, p->dn_vec_n, p->dn_vec_m, p->cusparse_handle, &p->buffer_size,
   &p->buffer);
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
