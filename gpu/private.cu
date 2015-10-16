#include "private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1
#define PRINT_INTERVAL 100

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)
#define NUM_SCALE_PASSES 1 /* additional passes don't help much */

#ifndef FLOAT
    #define CUBLAS(x) cublasD ## x
    #define CUSPARSE(x) cusparseD ## x
#else
    #define CUBLAS(x) cublasS ## x
    #define CUSPARSE(x) cusparseS ## x
#endif

static scs_int totCgIts;
static timer linsysTimer;
static scs_float totalSolveTime;

scs_int validateLinSys(const AMatrix * A) {
	scs_int i, rMax, Anz;
	if (!A->x || !A->i || !A->p) {
		scs_printf("data incompletely specified\n");
		return -1;
	}
	/* detects some errors in A col ptrs: */
	for (i = 0; i < A->n; ++i) {
		if (A->p[i] == A->p[i + 1]) {
			scs_printf("WARN: A->p (column pointers) not strictly increasing, column %li empty\n", (long) i);
		} else if (A->p[i] > A->p[i + 1]) {
			scs_printf("ERROR: A->p (column pointers) decreasing\n");
			return -1;
		}
	}
	Anz = A->p[A->n];
	if (((scs_float) Anz / A->m > A->n) || (Anz <= 0)) {
		scs_printf("Anz (nonzeros in A) = %li, outside of valid range\n", (long) Anz);
		return -1;
	}
	rMax = 0;
	for (i = 0; i < Anz; ++i) {
		if (A->i[i] > rMax)
			rMax = A->i[i];
	}
	if (rMax > A->m - 1) {
		scs_printf("number of rows in A inconsistent with input dimension\n");
		return -1;
	}
	return 0;
}

scs_int copyAMatrix(AMatrix ** dstp, const AMatrix * src) {
	scs_int Anz = src->p[src->n];
	AMatrix * A = (AMatrix *)scs_calloc(1, sizeof(AMatrix));
	if (!A) return 0;
	A->n = src->n;
	A->m = src->m;
	A->x = (scs_float *)scs_malloc(sizeof(scs_float) * Anz); /* A values, size: NNZ A */
	A->i = (scs_int *)scs_malloc(sizeof(scs_int) * Anz); /* A row index, size: NNZ A */
	A->p = (scs_int *)scs_malloc(sizeof(scs_int) * (src->n + 1)); /* A column pointer, size: n+1 */
	if (!A->x || !A->i || !A->p) return 0;
	memcpy(A->x, src->x, sizeof(scs_float) * Anz);
	memcpy(A->i, src->i, sizeof(scs_int) * Anz);
	memcpy(A->p, src->p, sizeof(scs_int) * (src->n + 1));
	*dstp = A;
	return 1;
}

void freeAMatrix(AMatrix * A) {
	if (A->x)
		scs_free(A->x);
	if (A->i)
		scs_free(A->i);
	if (A->p)
		scs_free(A->p);
	scs_free(A);
}

void printAMatrix(const AMatrix * A) {
	scs_int i, j;
	/* TODO: this is to prevent clogging stdout */
	if (A->p[A->n] < 2500) {
		scs_printf("\n");
		for (i = 0; i < A->n; ++i) {
			scs_printf("Col %li: ", (long) i);
			for (j = A->p[i]; j < A->p[i + 1]; j++) {
				scs_printf("A[%li,%li] = %4f, ", (long) A->i[j], (long) i, A->x[j]);
			}
			scs_printf("norm col = %4f\n", calcNorm(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
		}
		scs_printf("norm A = %4f\n", calcNorm(A->x, A->p[A->n]));
	}
}

void normalizeA(AMatrix * A, const Settings * stgs, const Cone * k, Scaling * scal) {
	scs_float * D = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
	scs_float * E = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
	scs_float * Dt = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
	scs_float * Et = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
	scs_float * nms = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
	scs_float minRowScale = MIN_SCALE * SQRTF((scs_float) A->n), maxRowScale = MAX_SCALE * SQRTF((scs_float) A->n);
	scs_float minColScale = MIN_SCALE * SQRTF((scs_float) A->m), maxColScale = MAX_SCALE * SQRTF((scs_float) A->m);
	scs_int i, j, l, count, delta, *boundaries, c1, c2;
	scs_float wrk, e;
	scs_int numBoundaries = getConeBoundaries(k, &boundaries);

#if EXTRAVERBOSE > 0
	timer normalizeTimer;
	tic(&normalizeTimer);
	scs_printf("normalizing A\n");
	printAMatrix(A);
#endif

	for (l = 0; l < NUM_SCALE_PASSES; ++l) {
		memset(D, 0, A->m * sizeof(scs_float));
		memset(E, 0, A->n * sizeof(scs_float));
		/* calculate row norms */
		for (i = 0; i < A->n; ++i) {
			c1 = A->p[i];
			c2 = A->p[i + 1];
			for (j = c1; j < c2; ++j) {
				wrk = A->x[j];
				D[A->i[j]] += wrk * wrk;
			}
		}
		for (i = 0; i < A->m; ++i) {
			D[i] = SQRTF(D[i]); /* just the norms */
		}

		/* mean of norms of rows across each cone  */
		count = boundaries[0];
		for (i = 1; i < numBoundaries; ++i) {
			wrk = 0;
			delta = boundaries[i];
			for (j = count; j < count + delta; ++j) {
				wrk += D[j];
			}
			wrk /= delta;
			for (j = count; j < count + delta; ++j) {
				D[j] = wrk;
			}
			count += delta;
		}

		for (i = 0; i < A->m; ++i) {
			if (D[i] < minRowScale)
				D[i] = 1;
			else if (D[i] > maxRowScale)
				D[i] = maxRowScale;
		}

		/* scale the rows with D */
		for (i = 0; i < A->n; ++i) {
			for (j = A->p[i]; j < A->p[i + 1]; ++j) {
				A->x[j] /= D[A->i[j]];
			}
		}
		/* calculate and scale by col norms, E */
		for (i = 0; i < A->n; ++i) {
			c1 = A->p[i + 1] - A->p[i];
			e = calcNorm(&(A->x[A->p[i]]), c1);
			if (e < minColScale)
				e = 1;
			else if (e > maxColScale)
				e = maxColScale;
			scaleArray(&(A->x[A->p[i]]), 1.0 / e, c1);
			E[i] = e;
		}

		for (i = 0; i < A->m; ++i) {
			Dt[i] = (l == 0) ? D[i] : Dt[i] * D[i];
		}
		for (i = 0; i < A->n; ++i) {
			Et[i] = (l == 0) ? E[i] : Et[i] * E[i];
		}
	}
	scs_free(boundaries);
	scs_free(D);
	scs_free(E);

	/* calculate mean of row norms of A */
	for (i = 0; i < A->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			wrk = A->x[j];
			nms[A->i[j]] += wrk * wrk;
		}
	}
	scal->meanNormRowA = 0.0;
	for (i = 0; i < A->m; ++i) {
		scal->meanNormRowA += SQRTF(nms[i]) / A->m;
	}
	scs_free(nms);

	/* calculate mean of col norms of A */
	scal->meanNormColA = 0.0;
	for (i = 0; i < A->n; ++i) {
		c1 = A->p[i + 1] - A->p[i];
		scal->meanNormColA += calcNorm(&(A->x[A->p[i]]), c1) / A->n;
	}

	/* scale up by d->SCALE if not equal to 1 */
	if (stgs->scale != 1) {
		scaleArray(A->x, stgs->scale, A->p[A->n]);
	}

	scal->D = Dt;
	scal->E = Et;

#if EXTRAVERBOSE > 0
	scs_printf("finished normalizing A, time: %1.2es\n", tocq(&normalizeTimer) / 1e3);
	printAMatrix(A);
#endif
}

void unNormalizeA(AMatrix * A, const Settings * stgs, const Scaling * scal) {
	scs_int i, j;
	scs_float * D = scal->D;
	scs_float * E = scal->E;
	for (i = 0; i < A->n; ++i) {
		scaleArray(&(A->x[A->p[i]]), E[i] / stgs->scale, A->p[i + 1] - A->p[i]);
	}
	for (i = 0; i < A->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			A->x[j] *= D[A->i[j]];
		}
	}
}

void accumByAtrans_gpu(const Priv * p, const scs_float *x, scs_float *y) {
    /* y += A'*x
       A in column compressed format
       parallelizes over columns (rows of A')
     */
    const scs_float onef = 1.0;
    AMatrix * Ag = p->Ag;
//    printf("ACCUM BY A TRANS GPU\n");
    CUSPARSE(csrmv)(p->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, Ag->n, Ag->m, p->Annz, &onef, p->descr, Ag->x, Ag->p, Ag->i, x, &onef, y); 
//    printf("FINISH ACCUM BY A TRANS GPU\n");
}

void accumByA_gpu(const Priv * p, const scs_float *x, scs_float *y) {
    /*y += A*x
      A in column compressed format
      this parallelizes over columns and uses
      pragma atomic to prevent concurrent writes to y
     */
    const scs_float onef = 1.0;
    AMatrix * Agt = p->Agt;
//    printf("ACCUM BY A GPU\n");
    CUSPARSE(csrmv)(p->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, Agt->n, Agt->m, p->Annz, &onef, p->descr, Agt->x, Agt->p, Agt->i, x, &onef, y);
//    printf("FINISH ACCUM BY A GPU\n");
}

void accumByAtrans(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
//    printf("ACCUM BY A TRANS\n");
    scs_float * v_m = p->tmp_m;
    scs_float * v_n = p->tmp_n;
    cudaMemcpy(v_m, x, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_n, y, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByAtrans_gpu(p, v_m, v_n);
    cudaMemcpy(y, v_n, A->n * sizeof(scs_float), cudaMemcpyDeviceToHost);
//    printf("FINISH ACCUM BY A TRANS\n");
}

void accumByA(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
//    printf("ACCUM BY A TRANS\n");
    scs_float * v_m = p->tmp_m;
    scs_float * v_n = p->tmp_n;
    cudaMemcpy(v_n, x, A->n * sizeof(scs_float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_m, y, A->m * sizeof(scs_float), cudaMemcpyHostToDevice);
    accumByA_gpu(p, v_n, v_m);
    cudaMemcpy(y, v_m, A->m * sizeof(scs_float), cudaMemcpyDeviceToHost);
//    printf("FINISH ACCUM BY A TRANS\n");
}

char * getLinSysMethod(const AMatrix * A, const Settings * s) {
	char * str = (char *)scs_malloc(sizeof(char) * 128);
	sprintf(str, "sparse-indirect GPU, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)", (long ) A->p[A->n], s->cg_rate);
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
    if(A->i)
        cudaFree(A->i);
    if(A->p)
        cudaFree(A->p);
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
		if (p->tmp_m)
			cudaFree(p->tmp_m);
        // XXX: do this with p later
		if (p->tmp_n)
			cudaFree(p->tmp_n);
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

/* solves (I+A'A)x = b, s warm start, solution stored in b */
/*y = (RHO_X * I + A'A)x */
static void matVec(const AMatrix * A, const Settings * s, Priv * p, const scs_float * x, scs_float * y) {
    /* x and y already loaded to GPU */
	scs_float * tmp_m = p->tmp_m;
	cudaMemset(tmp_m, 0, A->m * sizeof(scs_float));
	accumByA_gpu(p, x, tmp_m);
	cudaMemset(y, 0, A->n * sizeof(scs_float));
	accumByAtrans_gpu(p, tmp_m, y);
    CUBLAS(axpy)(p->cublasHandle, A->n, &(s->rho_x), x, 1, y, 1);
}

Priv * initPriv(const AMatrix * A, const Settings * stgs) {
	Priv * p = (Priv *)scs_calloc(1, sizeof(Priv));
    p->Annz = A->p[A->n];
    p->cublasHandle = 0;
    p->cusparseHandle = 0;
    p->descr = 0;

	totalSolveTime = 0;
	totCgIts = 0;

    /* Get handle to the CUBLAS context */
    p->cublasStatus = cublasCreate(&p->cublasHandle);
    //checkCudaErrors(p->cublasStatus);

    /* Get handle to the CUSPARSE context */
    p->cusparseStatus = cusparseCreate(&p->cusparseHandle);
    //checkCudaErrors(p->cusparseStatus);

    p->cusparseStatus = cusparseCreateMatDescr(&p->descr);
    //checkCudaErrors(p->cusparseStatus);

    cusparseSetMatType(p->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(p->descr, CUSPARSE_INDEX_BASE_ZERO);

    AMatrix * Ag = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Ag->n = A->n;
    Ag->m = A->m;
    p->Ag = Ag;

    AMatrix * Agt = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Agt->n = A->m;
    Agt->m = A->n;
    p->Agt = Agt;
   
    cudaMalloc((void **)&Ag->i, (A->p[A->n]) * sizeof(scs_int));
    cudaMalloc((void **)&Ag->p, (A->n + 1) * sizeof(scs_int));
    cudaMalloc((void **)&Ag->x, (A->p[A->n]) * sizeof(scs_float));

    cudaMalloc((void **)&p->p, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->r, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->Gp, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->bg, A->n * sizeof(scs_float));
    cudaMalloc((void **)&p->tmp_m, A->m * sizeof(scs_float)); /* intermediate result */
    cudaMalloc((void **)&p->tmp_n, A->n * sizeof(scs_float)); /* intermediate result */

    cudaMemcpy(Ag->i, A->i, (A->p[A->n]) * sizeof(scs_int), cudaMemcpyHostToDevice);
    cudaMemcpy(Ag->p, A->p, (A->n + 1) * sizeof(scs_int), cudaMemcpyHostToDevice);
    cudaMemcpy(Ag->x, A->x, (A->p[A->n]) * sizeof(scs_float), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&Agt->i, (A->p[A->n]) * sizeof(scs_int));
    cudaMalloc((void **)&Agt->p, (A->m + 1) * sizeof(scs_int));
    cudaMalloc((void **)&Agt->x, (A->p[A->n]) * sizeof(scs_float));

    /* transpose Ag into Agt */
    CUSPARSE(csr2csc)(p->cusparseHandle, A->n, A->m, A->p[A->n], Ag->x, Ag->p, Ag->i, Agt->x, Agt->i, Agt->p, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
    
    //CUDA_CHECK_ERR();
    
    if (!p->p || !p->r || !p->Gp || !p->tmp_m || !p->Ag || !p->Ag->i || !p->Ag->p || !p->Ag->x) {
		freePriv(p);
		return SCS_NULL;
	}
	return p;
}


static scs_int pcg(const AMatrix * A, const Settings * stgs, Priv * pr, const scs_float * s, scs_float * bh, scs_int max_its,
		scs_float tol) {
	scs_int i, n = A->n;
	scs_float alpha, nm_r, nm_r_old, ippGp, negAlpha, beta;
    scs_float onef = 1.0, negOnef = -1.0;
	scs_float *p = pr->p; /* cg direction */
	scs_float *Gp = pr->Gp; /* updated CG direction */
	scs_float *r = pr->r; /* cg residual */
    scs_float *b = pr->bg; /* rhs */
	cublasHandle_t cublasHandle = pr->cublasHandle;

    if (s == SCS_NULL) {
		cudaMemcpy(r, bh, n * sizeof(scs_float), cudaMemcpyHostToDevice);
  		cudaMemset(b, 0, n * sizeof(scs_float));
	} else {
        // bg (= b) contains s
        cudaMemcpy(b, s, n * sizeof(scs_float), cudaMemcpyHostToDevice);
        matVec(A, stgs, pr, b, r);

        // p contains bh temporarily
        cudaMemcpy(p, bh, n * sizeof(scs_float), cudaMemcpyHostToDevice);
        CUBLAS(axpy)(cublasHandle, n, &negOnef, p, 1, r, 1);
        CUBLAS(scal)(cublasHandle, n, &negOnef, r, 1);
    }
    
    CUBLAS(nrm2)(cublasHandle, n, r, 1, &nm_r);
//        CUBLAS(dot)(cublasHandle, n, r, 1, r, 1, &nm_r);
//        nm_r = SQRTF(nm_r);
    /* check to see if we need to run CG at all */
    if (nm_r < MIN(tol, 1e-18)) {
        scs_printf("early\n");
        return 0;
    }
    
    // put p in r, replacing temp mem
    cudaMemcpy(p, r, n * sizeof(scs_float), cudaMemcpyDeviceToDevice);

    for (i = 0; i < max_its; ++i) {
        matVec(A, stgs, pr, p, Gp);
        
        CUBLAS(dot)(cublasHandle, n, p, 1, Gp, 1, &ippGp);
        
        alpha = (nm_r * nm_r) / ippGp;
        negAlpha = -alpha;
        
        CUBLAS(axpy)(cublasHandle, n, &alpha, p, 1, b, 1);
        CUBLAS(axpy)(cublasHandle, n, &negAlpha, Gp, 1, r, 1);

        nm_r_old = nm_r;
//        CUBLAS(dot)(cublasHandle, n, r, 1, r, 1, &nm_r);
//        nm_r = SQRTF(nm_r);
        CUBLAS(nrm2)(cublasHandle, n, r, 1, &nm_r);
		if (nm_r < tol) {
#if EXTRAVERBOSE > 0
			scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, nm_r, (long) i+1);
#endif
            cudaMemcpy(bh, b, n * sizeof(scs_float), cudaMemcpyDeviceToHost);
			return i + 1;
		}
        beta = (nm_r * nm_r) / (nm_r_old * nm_r_old);
		CUBLAS(scal)(cublasHandle, n, &beta, p, 1);
        CUBLAS(axpy)(cublasHandle, n, &onef, r, 1, p, 1);
	}
    cudaMemcpy(bh, b, n * sizeof(scs_float), cudaMemcpyDeviceToHost);
    return i;
}

void accumByAtransO(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
	_accumByAtrans(A->n, A->x, A->i, A->p, x, y);
}

void accumByAO(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
	_accumByA(A->n, A->x, A->i, A->p, x, y);
}


void _accumByAtrans(scs_int n, scs_float * Ax, scs_int * Ai, scs_int * Ap, const scs_float *x, scs_float *y) {
    /* y += A'*x
       A in column compressed format
       parallelizes over columns (rows of A')
     */
    scs_int p, j;
    scs_int c1, c2;
    scs_float yj;
#if EXTRAVERBOSE > 0
    timer multByAtransTimer;
    tic(&multByAtransTimer);
#endif
#ifdef OPENMP
#pragma omp parallel for private(p,c1,c2,yj)
#endif
    for (j = 0; j < n; j++) {
        yj = y[j];
        c1 = Ap[j];
        c2 = Ap[j + 1];
        for (p = c1; p < c2; p++) {
            yj += Ax[p] * x[Ai[p]];
        }
        y[j] = yj;
    }
#if EXTRAVERBOSE > 0
    scs_printf("mult By A trans time: %1.2es\n", tocq(&multByAtransTimer) / 1e3);
#endif
}

void _accumByA(scs_int n, scs_float * Ax, scs_int * Ai, scs_int * Ap, const scs_float *x, scs_float *y) {
    /*y += A*x
      A in column compressed format
      this parallelizes over columns and uses
      pragma atomic to prevent concurrent writes to y
     */
    scs_int p, j;
    scs_int c1, c2;
    scs_float xj;
#if EXTRAVERBOSE > 0
    timer multByATimer;
    tic(&multByATimer);
#endif
    /*#pragma omp parallel for private(p,c1,c2,xj)  */
    for (j = 0; j < n; j++) {
        xj = x[j];
        c1 = Ap[j];
        c2 = Ap[j + 1];
        for (p = c1; p < c2; p++) {
            /*#pragma omp atomic */
            y[Ai[p]] += Ax[p] * xj;
        }
    }
#if EXTRAVERBOSE > 0
    scs_printf("mult By A time: %1.2es\n", tocq(&multByATimer) / 1e3);
#endif
}

scs_int solveLinSys(const AMatrix * A, const Settings * stgs, Priv * p, scs_float * b, const scs_float * s, scs_int iter) {
	scs_int cgIts;
	scs_float cgTol = calcNorm(b, A->n) * (iter < 0 ? CG_BEST_TOL : CG_MIN_TOL / POWF((scs_float) iter + 1, stgs->cg_rate));

	tic(&linsysTimer);
	/* solves Mx = b, for x but stores result in b */
	/* s contains warm-start (if available) */
	scs_float *tn, *tm;
    tn = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
    tm = (scs_float *)scs_malloc(A->m * sizeof(scs_float));

    scs_int n = A->n, m = A->m;
//    memcpy(tn, b, n * sizeof(scs_float));
//    memcpy(tm, &(b[A->n]), m * sizeof(scs_float));
    accumByAtrans(A, p, &(b[A->n]), b);
//    accumByAtransO(A, p, tm, tn);
//    printf("A trans err %4f\n", calcNormDiff(b, tn, n));

    /* solves (I+A'A)x = b, s warm start, solution stored in b */
	cgIts = pcg(A, stgs, p, s, b, A->n, MAX(cgTol, CG_BEST_TOL));
	scaleArray(&(b[A->n]), -1, A->m);

//    memcpy(tn, b, n * sizeof(scs_float));
//    memcpy(tm, &(b[A->n]), m * sizeof(scs_float));
    accumByA(A, p, b, &(b[A->n]));
//    accumByAO(A, p, tn, tm);
//    printf("A err %4f\n", calcNormDiff(&(b[A->n]), tm, m));

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

