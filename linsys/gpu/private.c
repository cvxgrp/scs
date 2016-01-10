#include "private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

#define CUDA_CHECK_ERR                                                         \
    do {                                                                       \
        cudaError_t err = cudaGetLastError();                                  \
        if (err != cudaSuccess) {                                              \
            printf("%s:%d:%s\n ERROR_CUDA: %s\n", __FILE__, __LINE__,          \
                   __func__, cudaGetErrorString(err));                         \
        }                                                                      \
    } while (0)

#ifndef EXTRAVERBOSE
#ifndef FLOAT
#define CUBLAS(x) cublasD##x
#define CUSPARSE(x) cusparseD##x
#else
#define CUBLAS(x) cublasS##x
#define CUSPARSE(x) cusparseS##x
#endif
#else
#ifndef FLOAT
#define CUBLAS(x)                                                              \
    CUDA_CHECK_ERR;                                                            \
    cublasD##x
#define CUSPARSE(x)                                                            \
    CUDA_CHECK_ERR;                                                            \
    cusparseD##x
#else
#define CUBLAS(x)                                                              \
    CUDA_CHECK_ERR;                                                            \
    cublasS##x
#define CUSPARSE(x)                                                            \
    CUDA_CHECK_ERR;                                                            \
    cusparseS##x
#endif
#endif

/*
 CUDA matrix routines only for CSR, not CSC matrices:
    CSC             CSR             GPU     Mult
    A (m x n)       A' (n x m)      Ag      accumByATransGpu
    A'(n x m)       A  (m x n)      Agt     accumByAGpu
*/

void accumByAtransGpu(const Priv *p, const scs_float *x, scs_float *y) {
    /* y += A'*x
       x and y MUST be on GPU already
    */
    const scs_float onef = 1.0;
    AMatrix *Ag = p->Ag;
    CUSPARSE(csrmv)(p->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, Ag->n,
                    Ag->m, p->Annz, &onef, p->descr, Ag->x, Ag->p, Ag->i, x,
                    &onef, y);
}

void accumByAGpu(const Priv *p, const scs_float *x, scs_float *y) {
    /* y += A*x
       x and y MUST be on GPU already
     */
    const scs_float onef = 1.0;
    AMatrix *Agt = p->Agt;
    CUSPARSE(csrmv)(p->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, Agt->n,
                    Agt->m, p->Annz, &onef, p->descr, Agt->x, Agt->p, Agt->i, x,
                    &onef, y);
}

/* do not use within pcg, reuses memory */
void accumByAtrans(const AMatrix *A, Priv *p, const scs_float *x,
                   scs_float *y) {
    scs_float *v_m = p->tmp_m;
    scs_float *v_n = p->r;
    cudaMemcpy(v_m, x, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_n, y, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByAtransGpu(p, v_m, v_n);
    cudaMemcpy(y, v_n, A->n * sizeof(scs_float), cudaMemcpyDeviceToHost);
}

/* do not use within pcg, reuses memory */
void accumByA(const AMatrix *A, Priv *p, const scs_float *x, scs_float *y) {
    scs_float *v_m = p->tmp_m;
    scs_float *v_n = p->r;
    cudaMemcpy(v_n, x, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_m, y, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByAGpu(p, v_n, v_m);
    cudaMemcpy(y, v_m, A->m * sizeof(scs_float), cudaMemcpyDeviceToHost);
}

char *getLinSysMethod(const AMatrix *A, const Settings *s) {
    char *str = (char *)scs_malloc(sizeof(char) * 128);
    sprintf(str, "sparse-indirect GPU, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)",
            (long)A->p[A->n], s->cg_rate);
    return str;
}

char *getLinSysSummary(Priv *p, const Info *info) {
    char *str = (char *)scs_malloc(sizeof(char) * 128);
    sprintf(str,
            "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
            (scs_float)p->totCgIts / (info->iter + 1),
            p->totalSolveTime / (info->iter + 1) / 1e3);
    p->totCgIts = 0;
    p->totalSolveTime = 0;
    return str;
}

void cudaFreeAMatrix(AMatrix *A) {
    if (A->x)
        cudaFree(A->x);
    if (A->i)
        cudaFree(A->i);
    if (A->p)
        cudaFree(A->p);
}

void freePriv(Priv *p) {
    if (p) {
        if (p->p)
            cudaFree(p->p);
        if (p->r)
            cudaFree(p->r);
        if (p->Gp)
            cudaFree(p->Gp);
        if (p->bg)
            cudaFree(p->bg);
        if (p->tmp_m)
            cudaFree(p->tmp_m);
        if (p->z)
            cudaFree(p->z);
        if (p->M)
            cudaFree(p->M);
        if (p->Ag) {
            cudaFreeAMatrix(p->Ag);
            scs_free(p->Ag);
        }
        if (p->Agt) {
            cudaFreeAMatrix(p->Agt);
            scs_free(p->Agt);
        }
        cusparseDestroy(p->cusparseHandle);
        cublasDestroy(p->cublasHandle);
        cudaDeviceReset();
        scs_free(p);
    }
}

/*y = (RHO_X * I + A'A)x */
static void matVec(const AMatrix *A, const Settings *s, Priv *p,
                   const scs_float *x, scs_float *y) {
    /* x and y MUST already be loaded to GPU */
    scs_float *tmp_m = p->tmp_m; /* temp memory */
    cudaMemset(tmp_m, 0, A->m * sizeof(scs_float));
    accumByAGpu(p, x, tmp_m);
    cudaMemset(y, 0, A->n * sizeof(scs_float));
    accumByAtransGpu(p, tmp_m, y);
    CUBLAS(axpy)(p->cublasHandle, A->n, &(s->rho_x), x, 1, y, 1);
}

/* M = inv ( diag ( RHO_X * I + A'A ) ) */
void getPreconditioner(const AMatrix *A, const Settings *stgs, Priv *p) {
    scs_int i;
    scs_float *M = (scs_float *)scs_malloc(A->n * sizeof(scs_float));

#if EXTRAVERBOSE > 0
    scs_printf("getting pre-conditioner\n");
#endif

    for (i = 0; i < A->n; ++i) {
        M[i] = 1 / (stgs->rho_x +
                    calcNormSq(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
        /* M[i] = 1; */
    }
    cudaMemcpy(p->M, M, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    scs_free(M);

#if EXTRAVERBOSE > 0
    scs_printf("finished getting pre-conditioner\n");
#endif
}

Priv *initPriv(const AMatrix *A, const Settings *stgs) {
    cudaError_t err;
    Priv *p = (Priv *)scs_calloc(1, sizeof(Priv));
    p->Annz = A->p[A->n];
    p->cublasHandle = 0;
    p->cusparseHandle = 0;
    p->descr = 0;

    p->totalSolveTime = 0;
    p->totCgIts = 0;

    /* Get handle to the CUBLAS context */
    cublasCreate(&p->cublasHandle);

    /* Get handle to the CUSPARSE context */
    cusparseCreate(&p->cusparseHandle);

    /* Matrix description */
    cusparseCreateMatDescr(&p->descr);
    cusparseSetMatType(p->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(p->descr, CUSPARSE_INDEX_BASE_ZERO);

    AMatrix *Ag = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Ag->n = A->n;
    Ag->m = A->m;
    p->Ag = Ag;

    AMatrix *Agt = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Agt->n = A->m;
    Agt->m = A->n;
    p->Agt = Agt;

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
    cudaMemcpy(Ag->p, A->p, (A->n + 1) * sizeof(scs_int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(Ag->x, A->x, (A->p[A->n]) * sizeof(scs_float),
               cudaMemcpyHostToDevice);

    cudaMalloc((void **)&Agt->i, (A->p[A->n]) * sizeof(scs_int));
    cudaMalloc((void **)&Agt->p, (A->m + 1) * sizeof(scs_int));
    cudaMalloc((void **)&Agt->x, (A->p[A->n]) * sizeof(scs_float));

    getPreconditioner(A, stgs, p);

    /* transpose Ag into Agt for faster multiplies */
    /* TODO: memory intensive, could perform transpose in CPU and copy to GPU */
    CUSPARSE(csr2csc)(p->cusparseHandle, A->n, A->m, A->p[A->n], Ag->x, Ag->p,
                      Ag->i, Agt->x, Agt->i, Agt->p, CUSPARSE_ACTION_NUMERIC,
                      CUSPARSE_INDEX_BASE_ZERO);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("%s:%d:%s\nERROR_CUDA: %s\n", __FILE__, __LINE__, __func__,
               cudaGetErrorString(err));
        freePriv(p);
        return SCS_NULL;
    }
    return p;
}

static void applyPreConditioner(cublasHandle_t cublasHandle, scs_float *M,
                                scs_float *z, scs_float *r, scs_int n) {
    cudaMemcpy(z, r, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    CUBLAS(tbmv)(cublasHandle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
                 CUBLAS_DIAG_NON_UNIT, n, 0, M, 1, z, 1);
}

/* solves (I+A'A)x = b, s warm start, solution stored in bg (on GPU) */
static scs_int pcg(const AMatrix *A, const Settings *stgs, Priv *pr,
                   const scs_float *s, scs_float *bg, scs_int max_its,
                   scs_float tol) {
    scs_int i, n = A->n;
    scs_float alpha, nrm_r, pGp, negAlpha, beta, ipzr, ipzrOld;
    scs_float onef = 1.0, negOnef = -1.0;
    scs_float *p = pr->p;   /* cg direction */
    scs_float *Gp = pr->Gp; /* updated CG direction */
    scs_float *r = pr->r;   /* cg residual */
    scs_float *z = pr->z;   /* preconditioned */
    scs_float *M = pr->M;   /* preconditioner */
    cublasHandle_t cublasHandle = pr->cublasHandle;

    if (s == SCS_NULL) {
        cudaMemcpy(r, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
        cudaMemset(bg, 0, n * sizeof(scs_float));
    } else {
        /* p contains bg temporarily */
        cudaMemcpy(p, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
        /* bg contains s */
        cudaMemcpy(bg, s, n * sizeof(scs_float), cudaMemcpyHostToDevice);
        matVec(A, stgs, pr, bg, r);
        CUBLAS(axpy)(cublasHandle, n, &negOnef, p, 1, r, 1);
        CUBLAS(scal)(cublasHandle, n, &negOnef, r, 1);
    }

    /* for some reason nrm2 is VERY slow */
    /* CUBLAS(nrm2)(cublasHandle, n, r, 1, &nrm_r); */
    CUBLAS(dot)(cublasHandle, n, r, 1, r, 1, &nrm_r);
    nrm_r = SQRTF(nrm_r);
    /* check to see if we need to run CG at all */
    if (nrm_r < MIN(tol, 1e-18)) {
        return 0;
    }

    applyPreConditioner(cublasHandle, M, z, r, n);
    CUBLAS(dot)(cublasHandle, n, r, 1, z, 1, &ipzr);
    /* put z in p, replacing temp mem */
    cudaMemcpy(p, z, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);

    for (i = 0; i < max_its; ++i) {
        matVec(A, stgs, pr, p, Gp);

        CUBLAS(dot)(cublasHandle, n, p, 1, Gp, 1, &pGp);

        alpha = ipzr / pGp;
        negAlpha = -alpha;

        CUBLAS(axpy)(cublasHandle, n, &alpha, p, 1, bg, 1);
        CUBLAS(axpy)(cublasHandle, n, &negAlpha, Gp, 1, r, 1);

        /* for some reason nrm2 is VERY slow */
        /* CUBLAS(nrm2)(cublasHandle, n, r, 1, &nrm_r); */
        CUBLAS(dot)(cublasHandle, n, r, 1, r, 1, &nrm_r);
        nrm_r = SQRTF(nrm_r);
        if (nrm_r < tol) {
            i++;
            break;
        }
        ipzrOld = ipzr;
        applyPreConditioner(cublasHandle, M, z, r, n);
        CUBLAS(dot)(cublasHandle, n, r, 1, z, 1, &ipzr);

        beta = ipzr / ipzrOld;
        CUBLAS(scal)(cublasHandle, n, &beta, p, 1);
        CUBLAS(axpy)(cublasHandle, n, &onef, z, 1, p, 1);
    }
#if EXTRAVERBOSE > 0
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, nrm_r, (long)i + 1);
#endif
    return i;
}

#ifdef TEST_GPU_MAT_MUL
void accumByAtransHost(const AMatrix *A, Priv *p, const scs_float *x,
                       scs_float *y) {
    _accumByAtrans(A->n, A->x, A->i, A->p, x, y);
}

void accumByAHost(const AMatrix *A, Priv *p, const scs_float *x, scs_float *y) {
    _accumByA(A->n, A->x, A->i, A->p, x, y);
}

void testGpuMatMul(const AMatrix *A, Priv *p, scs_float *b) {
    /* test to see if matrix multiplication codes agree */
    scs_float t[A->n + A->m], u[A->n + A->m], *bg;
    cudaMalloc((void **)&bg, (A->n + A->m) * sizeof(scs_float));

    cudaMemcpy(bg, b, (A->n + A->m) * sizeof(scs_float),
               cudaMemcpyHostToDevice);
    memcpy(t, b, (A->n + A->m) * sizeof(scs_float));

    accumByAtransGpu(p, &(bg[A->n]), bg);
    accumByAtransHost(A, p, &(t[A->n]), t);
    cudaMemcpy(u, bg, (A->n + A->m) * sizeof(scs_float),
               cudaMemcpyDeviceToHost);
    printf("A trans multiplication err %2.e\n", calcNormDiff(u, t, A->n));

    accumByAGpu(p, bg, &(bg[A->n]));
    accumByAHost(A, p, t, &(t[A->n]));
    cudaMemcpy(u, bg, (A->n + A->m) * sizeof(scs_float),
               cudaMemcpyDeviceToHost);
    printf("A multiplcation err %2.e\n",
           calcNormDiff(&(u[A->n]), &(t[A->n]), A->m));
    cudaFree(bg);
}
#endif

scs_int solveLinSys(const AMatrix *A, const Settings *stgs, Priv *p,
                    scs_float *b, const scs_float *s, scs_int iter) {
    scs_int cgIts;
    timer linsysTimer;
    scs_float *bg = p->bg;
    scs_float negOnef = -1.0;
    scs_float cgTol =
        calcNorm(b, A->n) *
        (iter < 0 ? CG_BEST_TOL
                  : CG_MIN_TOL / POWF((scs_float)iter + 1, stgs->cg_rate));

    tic(&linsysTimer);
/* solves Mx = b, for x but stores result in b */
/* s contains warm-start (if available) */

#ifdef TEST_GPU_MAT_MUL
    testGpuMatMul(A, p, b);
#endif

    /* all on GPU */
    cudaMemcpy(bg, b, (A->n + A->m) * sizeof(scs_float),
               cudaMemcpyHostToDevice);
    accumByAtransGpu(p, &(bg[A->n]), bg);
    /* solves (I+A'A)x = b, s warm start, solution stored in b */
    cgIts = pcg(A, stgs, p, s, bg, A->n, MAX(cgTol, CG_BEST_TOL));
    CUBLAS(scal)(p->cublasHandle, A->m, &negOnef, &(bg[A->n]), 1);
    accumByAGpu(p, bg, &(bg[A->n]));
    cudaMemcpy(b, bg, (A->n + A->m) * sizeof(scs_float),
               cudaMemcpyDeviceToHost);

    if (iter >= 0) {
        p->totCgIts += cgIts;
    }

    p->totalSolveTime += tocq(&linsysTimer);
#if EXTRAVERBOSE > 0
    scs_printf("linsys solve time: %1.2es\n", tocq(&linsysTimer) / 1e3);
#endif
    return 0;
}

#ifdef __cplusplus
}
#endif
