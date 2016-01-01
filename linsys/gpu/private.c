#include "private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

#define CUDA_CHECK_ERR \
  do { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
      printf("%s:%d:%s\n ERROR_CUDA: %s\n", __FILE__, __LINE__, __func__, \
             cudaGetErrorString(err)); \
    } \
  } while (0)

#ifndef EXTRAVERBOSE
#ifndef FLOAT
    #define CUBLAS(x) cublasD ## x
#else
    #define CUBLAS(x) cublasS ## x
#endif
#else
#ifndef FLOAT
    #define CUBLAS(x) CUDA_CHECK_ERR; cublasD ## x
#else
    #define CUBLAS(x) CUDA_CHECK_ERR; cublasS ## x
#endif
#endif

static scs_int totCgIts;
static timer linsysTimer;
static scs_float totalSolveTime;

void accumByAtransGpu(const Priv * p, const scs_float *x, scs_float *y) {
    /* y += A'*x
       x and y MUST be on GPU already
    */
    const scs_float onef = 1.0;
    AMatrix * Ag = p->Ag;
    CUBLAS(gemv)(p->cublasHandle, CUBLAS_OP_T, Ag->m, Ag->n, &onef, Ag->x, Ag->m, x, 1, &onef, y, 1);
}

void accumByAGpu(const Priv * p, const scs_float *x, scs_float *y) {
    /* y += A*x
       x and y MUST be on GPU already
     */
    const scs_float onef = 1.0;
    AMatrix * Ag = p->Ag;
    CUBLAS(gemv)(p->cublasHandle, CUBLAS_OP_N, Ag->m, Ag->n, &onef, Ag->x, Ag->m, x, 1, &onef, y, 1);
}

/* do not use within pcg, reuses memory */
void accumByAtrans(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
    scs_float * v_m = p->bg;
    scs_float * v_n = p->r;
    cudaMemcpy(v_m, x, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_n, y, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByAtransGpu(p, v_m, v_n);
    cudaMemcpy(y, v_n, A->n * sizeof(scs_float), cudaMemcpyDeviceToHost);
}

/* do not use within pcg, reuses memory */
void accumByA(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
    scs_float * v_m = p->bg;
    scs_float * v_n = p->r;
    cudaMemcpy(v_n, x, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_m, y, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByAGpu(p, v_n, v_m);
    cudaMemcpy(y, v_m, A->m * sizeof(scs_float), cudaMemcpyDeviceToHost);
}

char * getLinSysMethod(const AMatrix * A, const Settings * stgs) {
    char * str = (char *)scs_malloc(sizeof(char) * 128);
    sprintf(str, "dense-indirect GPU, CG tol ~ 1/iter^(%2.2f)", stgs->cg_rate);
    return str;
}

char * getLinSysSummary(Priv * p, const Info * info) {
	char * str = (char *)scs_malloc(sizeof(char) * 128);
	sprintf(str, "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
			(scs_float ) totCgIts / (info->iter + 1), totalSolveTime / (info->iter + 1) / 1e3);
	totCgIts = 0;
	totalSolveTime = 0;
	return str;
}

void cudaFreeAMatrix(AMatrix * A) {
    if(A->x)
        cudaFree(A->x);
}

void freePriv(Priv * p) {
    if (p) {
		if (p->p)
			cudaFree(p->p);
		if (p->r)
			cudaFree(p->r);
		if (p->Gp)
			cudaFree(p->Gp);
		if (p->bg)
			cudaFree(p->bg);
		if (p->z)
			cudaFree(p->z);
		if (p->M)
			cudaFree(p->M);
		if (p->Ag) {
            cudaFreeAMatrix(p->Ag);
			scs_free(p->Ag);
		}
		if (p->G) {
            cudaFreeAMatrix(p->G);
			scs_free(p->G);
		}
        cublasDestroy(p->cublasHandle);
        cudaDeviceReset();
        scs_free(p);
	}
}

/* store inverse preconditioner */
void getPreconditioner(Priv *p) {
    scs_float onef = 1.0;
    scs_int i, n = p->G->n;
    scs_float * tmp = scs_malloc(n * sizeof(scs_float));
    
    cudaMemset(p->M, 0, n * sizeof(scs_float));
    CUBLAS(axpy)(p->cublasHandle, n, &onef, p->G->x, n + 1, p->M, 1);

    /* annoying hack to invert M */
    cudaMemcpy(tmp, p->M, n * sizeof(scs_float), cudaMemcpyDeviceToHost);
    for(i = 0; i < n; i++) {
        tmp[i] = 1 / tmp[i];
    }
    cudaMemcpy(p->M, tmp, n * sizeof(scs_float), cudaMemcpyHostToDevice);
    scs_free(tmp);
}

Priv * initPriv(const AMatrix * A, const Settings * stgs) {
    Priv * p = (Priv *)scs_calloc(1, sizeof(Priv));
    scs_float *Gtemp, onef = 1.0;
    cudaError_t err;
    scs_int j, numElementsG;

    p->cublasHandle = 0;
	totalSolveTime = 0;
	totCgIts = 0;

    /* Get handle to the CUBLAS context */
    cublasCreate(&p->cublasHandle);

    AMatrix * Ag = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Ag->n = A->n;
    Ag->m = A->m;
    p->Ag = Ag;
    AMatrix * G = (AMatrix *)scs_malloc(sizeof(AMatrix));
	G->n = A->n;
	G->m = A->n;
	p->G = G;
	numElementsG = A->n * A->n;

    cudaMalloc((void **)&Ag->x, A->n * A->m * sizeof(scs_float));
    cudaMalloc((void **)&G->x, numElementsG * sizeof(scs_float));

    cudaMalloc((void **)&p->p, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->r, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->Gp, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->bg, (A->n + A->m) * sizeof(scs_float));
    cudaMalloc((void **)&p->z, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->M, A->n * sizeof(scs_float));

    cudaMemcpy(Ag->x, A->x, A->n * A->m * sizeof(scs_float), cudaMemcpyHostToDevice);

    /* create G in CPU mem first (just rho_x * I part), then send to GPU */
	Gtemp = (scs_float *)scs_calloc(A->n * A->n, sizeof(scs_float));
    for (j = 0; j < A->n; j++) {
        Gtemp[j * A->n + j] = stgs->rho_x;
    }
    cudaMemcpy(G->x, Gtemp, A->n * A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    scs_free(Gtemp);

	CUBLAS(syrk)(p->cublasHandle, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, A->n, A->m, &onef, p->Ag->x, A->m, &onef, G->x, A->n);
    getPreconditioner(p);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("%s:%d:%s\nERROR_CUDA: %s\n", __FILE__, __LINE__, __func__, cudaGetErrorString(err));
		freePriv(p);
		return SCS_NULL;
	}
	return p;
}

/* use actual preconditioner (not inverse) */
static void applyPreConditioner(cublasHandle_t cublasHandle, scs_float * M, scs_float * z, scs_float * r, scs_int n) {
    cudaMemcpy(z, r, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    CUBLAS(tbmv)(cublasHandle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, 0, M, 1, z, 1);
}

/* solves (I+A'A)x = b, s warm start, solution stored in bg (on GPU) */
static scs_int pcg(const AMatrix * A, const Settings * stgs, Priv * pr, const scs_float * s, scs_float * bg, scs_int max_its,
        scs_float tol) {
    scs_int i, n = A->n;
    scs_float alpha, nrm_r, pGp, negAlpha, beta, ipzr, ipzrOld;
    scs_float onef = 1.0, negOnef = -1.0, zerof = 0.0;
    scs_float *p = pr->p; /* cg direction */
    scs_float *Gp = pr->Gp; /* updated CG direction */
    scs_float *r = pr->r; /* cg residual */
    scs_float *z = pr->z; /* preconditioned */
    scs_float *M = pr->M; /* preconditioner */
    cublasHandle_t cublasHandle = pr->cublasHandle;

	cudaMemcpy(r, bg, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);
    if (s == SCS_NULL) {
        cudaMemset(bg, 0, n * sizeof(scs_float));
    } else {
        /* bg contains s */
        cudaMemcpy(bg, s, n * sizeof(scs_float), cudaMemcpyHostToDevice);
		CUBLAS(symv)(pr->cublasHandle, CUBLAS_FILL_MODE_UPPER, n, &negOnef, pr->G->x, n, bg, 1, &onef, r, 1);
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
    	CUBLAS(symv)(pr->cublasHandle, CUBLAS_FILL_MODE_UPPER, n, &onef, pr->G->x, n, p, 1, &zerof, Gp, 1);

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
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, nrm_r, (long) i+1);
#endif
    return i;
}

#ifdef TEST_GPU_MAT_MUL
void accumByAtransHost(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
    _accumByAtrans(A->n, A->x, A->i, A->p, x, y);
}

void accumByAHost(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
    _accumByA(A->n, A->x, A->i, A->p, x, y);
}

void testGpuMatMul(const AMatrix * A, Priv * p, scs_float * b) {
    /* test to see if matrix multiplication codes agree */
    scs_float t[A->n + A->m], u[A->n + A->m], *bg;
    cudaMalloc((void **)&bg, (A->n + A->m) * sizeof(scs_float));

    cudaMemcpy(bg, b, (A->n + A->m) * sizeof(scs_float), cudaMemcpyHostToDevice);
    memcpy(t, b, (A->n + A->m) * sizeof(scs_float));

    accumByAtransGpu(p, &(bg[A->n]), bg);
    accumByAtransHost(A, p, &(t[A->n]), t);
    cudaMemcpy(u, bg, (A->n + A->m) * sizeof(scs_float), cudaMemcpyDeviceToHost);
    printf("A trans multiplication err %2.e\n", calcNormDiff(u, t, A->n));

    accumByAGpu(p, bg, &(bg[A->n]));
    accumByAHost(A, p, t, &(t[A->n]));
    cudaMemcpy(u, bg, (A->n + A->m) * sizeof(scs_float), cudaMemcpyDeviceToHost);
    printf("A multiplcation err %2.e\n", calcNormDiff(&(u[A->n]), &(t[A->n]), A->m));
    cudaFree(bg);
}
#endif

scs_int solveLinSys(const AMatrix * A, const Settings * stgs, Priv * p, scs_float * b, const scs_float * s, scs_int iter) {
    scs_int cgIts;
    scs_float * bg = p->bg;
    scs_float negOnef = -1.0;
    scs_float cgTol = calcNorm(b, A->n) * (iter < 0 ? CG_BEST_TOL : CG_MIN_TOL / POWF((scs_float) iter + 1, stgs->cg_rate));

    tic(&linsysTimer);
    /* solves Mx = b, for x but stores result in b */
    /* s contains warm-start (if available) */

#ifdef TEST_GPU_MAT_MUL
    testGpuMatMul(A, p, b);
#endif

    /* all on GPU */
    cudaMemcpy(bg, b, (A->n + A->m) * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByAtransGpu(p, &(bg[A->n]), bg);
    /* solves (I+A'A)x = b, s warm start, solution stored in b */
    cgIts = pcg(A, stgs, p, s, bg, A->n, MAX(cgTol, CG_BEST_TOL));
    CUBLAS(scal)(p->cublasHandle, A->m, &negOnef, &(bg[A->n]), 1);
    accumByAGpu(p, bg, &(bg[A->n]));
    cudaMemcpy(b, bg, (A->n + A->m) * sizeof(scs_float), cudaMemcpyDeviceToHost);

    if (iter >= 0) {
        totCgIts += cgIts;
    }

    totalSolveTime += tocq(&linsysTimer);
#if EXTRAVERBOSE > 0
    scs_printf("linsys solve time: %1.2es\n", tocq(&linsysTimer) / 1e3);
#endif
    return 0;
}

#ifdef __cplusplus
}
#endif

