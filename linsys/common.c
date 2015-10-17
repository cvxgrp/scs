#include "common.h"

scs_int validateLinSys(const AMatrix * A) {
	return (A && A->x);
}

scs_int copyAMatrix(AMatrix ** dstp, const AMatrix * src) {
   AMatrix * A = scs_calloc(1, sizeof(AMatrix));
   if (!A) return 0;
   A->n = src->n;
   A->m = src->m;
   A->x = scs_malloc(sizeof(scs_float) * src->n * src->m);
   if (!A->x) return 0;
   memcpy(A->x, src->x, sizeof(scs_float) * src->n * src->m);
   *dstp = A;
   return 1; 
}

void freeAMatrix(AMatrix * A) {
    if (A->x)
        scs_free(A->x);
    scs_free(A);
}

void printAMatrix(AMatrix * A) {
    scs_int i, j;
    /* TODO: this is to prevent clogging stdout */
    if (A->n * A->m < 2500) {
        scs_printf("\n");
        for (i = 0; i < A->n; ++i) {
            scs_printf("Col %li: ", (long) i);
            for (j = 0; j < A->m; ++j) {
                scs_printf("A[%li,%li] = %4f, ", (long) j, (long) i, A->x[i * A->m + j]);
            }
            scs_printf("norm col = %4f\n", calcNorm(&(A->x[i * A->m]), A->m));
        }
        scs_printf("norm A = %4f\n", calcNorm(A->x, A->n * A->m));
    }
}

void normalizeA(AMatrix * A, const Settings * stgs, const Cone * k, Scaling * scal) {
	scs_int i, j, count, delta;
	scs_float wrk;
    scs_float * D = scs_calloc(A->m, sizeof(scs_float));
    scs_float * E = scs_calloc(A->n, sizeof(scs_float));
    scs_float minRowScale = MIN_SCALE * SQRTF(A->n), maxRowScale = MAX_SCALE * SQRTF(A->n);
    scs_float minColScale = MIN_SCALE * SQRTF(A->m), maxColScale = MAX_SCALE * SQRTF(A->m);
    scs_int *boundaries, numBoundaries = getConeBoundaries(k, &boundaries);
    blasint n = (blasint) A->n, m = (blasint) A->m, one = 1, nm = n * m;

#ifdef EXTRAVERBOSE
    timer normalizeTimer;
    tic(&normalizeTimer);
    scs_printf("normalizing A\n");
    printAMatrix(A);
#endif

	/* calculate row norms */
	for (i = 0; i < A->m; ++i) {
		D[i] = BLAS(nrm2)(&n, &(A->x[i]), &m);
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
	scs_free(boundaries);

    for (i = 0; i < A->m; ++i) {
        if (D[i] < minRowScale)
            D[i] = 1;
        else if (D[i] > maxRowScale)
            D[i] = maxRowScale;
    }

    /* scale the rows with D */
	for (i = 0; i < A->m; ++i) {
		wrk = 1.0 / D[i];
		BLAS(scal)(&n, &wrk, &(A->x[i]), &m);
	}

	/* calculate and scale by col norms, E */
	for (i = 0; i < A->n; ++i) {
        E[i] = BLAS(nrm2)(&m, &(A->x[i * A->m]), &one);
        if (E[i] < minColScale)
            E[i] = 1;
        else if (E[i] > maxColScale)
            E[i] = maxColScale;
        wrk = 1.0 / E[i];
        BLAS(scal)(&m, &wrk, &(A->x[i * A->m]), &one);
    }

	scal->meanNormRowA = 0.0;
	for (i = 0; i < A->m; ++i) {
		scal->meanNormRowA += BLAS(nrm2)(&n, &(A->x[i]), &m) / A->m;
	}

    scal->meanNormColA = 0.0;
	for (i = 0; i < A->n; ++i) {
		scal->meanNormColA += BLAS(nrm2)(&m, &(A->x[i * A->m]), &one) / A->n;
	}

	if (stgs->scale != 1) {
		BLAS(scal)(&nm, &stgs->scale, A->x, &one);
	}

	scal->D = D;
	scal->E = E;

#ifdef EXTRAVERBOSE
    scs_printf("finished normalizing A, time: %1.2es\n", tocq(&normalizeTimer) / 1e3);
    printAMatrix(A);
#endif
}

void unNormalizeA(AMatrix * A, const Settings * stgs, const Scaling * scal) {
	scs_int i;
	scs_float * D = scal->D;
	scs_float * E = scal->E;
	scs_float invScale = 1.0 / stgs->scale;
	blasint n = (blasint) A->n, m = (blasint) A->m, one = 1, nm = n * m;
	for (i = 0; i < A->n; ++i) {
		BLAS(scal)(&m, &(E[i]), &(A->x[i * A->m]), &one);
	}
	for (i = 0; i < A->m; ++i) {
		BLAS(scal)(&n, &(D[i]), &(A->x[i]), &m);
	}
	BLAS(scal)(&nm, &invScale, A->x, &one);
}

void _accumByAtrans(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
    blasint one = 1, n = (blasint) A->n, m = (blasint) A->m;
    scs_float onef = 1.0;
    BLAS(gemv)("Transpose", &m, &n, &onef, A->x, &m, x, &one, &onef, y, &one);
}

void _accumByA(const AMatrix * A, Priv * p, const scs_float *x, scs_float *y) {
    blasint one = 1, n = (blasint) A->n, m = (blasint) A->m;
    scs_float onef = 1.0;
    BLAS(gemv)("NoTranspose", &m, &n, &onef, A->x, &m, x, &one, &onef, y, &one);
}

