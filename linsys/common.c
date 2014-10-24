#include "common.h"

scs_int validateLinSys(Data *d) {
	return (d->A && d->A->x);
}

void freeAMatrix(AMatrix * A) {
    if (A->x)
        scs_free(A->x);
}

void printAMatrix(Data * d) {
    scs_int i, j;
    AMatrix * A = d->A;
    /* TODO: this is to prevent clogging stdout */
    if (d->n * d->m < 2500) {
        scs_printf("\n");
        for (i = 0; i < d->n; ++i) {
            scs_printf("Col %li: ", (long) i);
            for (j = 0; j < d->m; ++j) {
                scs_printf("A[%li,%li] = %4f, ", (long) j, (long) i, A->x[i * d->m + j]);
            }
            scs_printf("norm col = %4f\n", calcNorm(&(A->x[i * d->m]), d->m));
        }
        scs_printf("norm A = %4f\n", calcNorm(A->x, d->n * d->m));
    }
}

void normalizeA(Data * d, Work * w, Cone * k) {
	scs_int i, j, count, delta;
	scs_float wrk;
    scs_float * D = scs_calloc(d->m, sizeof(scs_float));
    scs_float * E = scs_calloc(d->n, sizeof(scs_float));
    scs_float minRowScale = MIN_SCALE * SQRTF(d->n), maxRowScale = MAX_SCALE * SQRTF(d->n);
    scs_float minColScale = MIN_SCALE * SQRTF(d->m), maxColScale = MAX_SCALE * SQRTF(d->m);
    scs_int *boundaries, numBoundaries = getConeBoundaries(k, &boundaries);
    blasint n = (blasint) d->n, m = (blasint) d->m, one = 1, nm = n * m;

#ifdef EXTRAVERBOSE
    timer normalizeTimer;
    tic(&normalizeTimer);
    scs_printf("normalizing A\n");
    printAMatrix(d);
#endif

	/* calculate row norms */
	for (i = 0; i < d->m; ++i) {
		D[i] = BLAS(nrm2)(&n, &(d->A->x[i]), &m);
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

    for (i = 0; i < d->m; ++i) {
        if (D[i] < minRowScale)
            D[i] = 1;
        else if (D[i] > maxRowScale)
            D[i] = maxRowScale;
    }

    /* scale the rows with D */
	for (i = 0; i < d->m; ++i) {
		wrk = 1.0 / D[i];
		BLAS(scal)(&n, &wrk, &(d->A->x[i]), &m);
	}

	/* calculate and scale by col norms, E */
	for (i = 0; i < d->n; ++i) {
        E[i] = BLAS(nrm2)(&m, &(d->A->x[i * d->m]), &one);
        if (E[i] < minColScale)
            E[i] = 1;
        else if (E[i] > maxColScale)
            E[i] = maxColScale;
        wrk = 1.0 / E[i];
        BLAS(scal)(&m, &wrk, &(d->A->x[i * d->m]), &one);
    }

	w->meanNormRowA = 0.0;
	for (i = 0; i < d->m; ++i) {
		w->meanNormRowA += BLAS(nrm2)(&n, &(d->A->x[i]), &m) / d->m;
	}

    w->meanNormColA = 0.0;
	for (i = 0; i < d->n; ++i) {
		w->meanNormColA += BLAS(nrm2)(&m, &(d->A->x[i * d->m]), &one) / d->n;
	}

	if (d->scale != 1) {
		BLAS(scal)(&nm, &d->scale, d->A->x, &one);
	}

	w->D = D;
	w->E = E;

#ifdef EXTRAVERBOSE
    scs_printf("finished normalizing A, time: %1.2es\n", tocq(&normalizeTimer) / 1e3);
    printAMatrix(d);
#endif
}

void unNormalizeA(Data *d, Work * w) {
	scs_int i;
	scs_float * D = w->D;
	scs_float * E = w->E;
	scs_float invScale = 1.0 / d->scale;
	blasint n = (blasint) d->n, m = (blasint) d->m, one = 1, nm = n * m;
	for (i = 0; i < d->n; ++i) {
		BLAS(scal)(&m, &(E[i]), &(d->A->x[i * d->m]), &one);
	}
	for (i = 0; i < d->m; ++i) {
		BLAS(scal)(&n, &(D[i]), &(d->A->x[i]), &m);
	}
	BLAS(scal)(&nm, &invScale, d->A->x, &one);
}

void accumByAtrans(Data * d, Priv * p, const scs_float * x, scs_float * y) {
	scs_float * A = d->A->x;
	blasint one = 1, n = (blasint) d->n, m = (blasint) d->m;
	scs_float onef = 1.0;
	BLAS(gemv)("Transpose", &m, &n, &onef, A, &m, x, &one, &onef, y, &one);
}

void accumByA(Data * d, Priv * p, const scs_float * x, scs_float * y) {
	scs_float * A = d->A->x;
	blasint one = 1, n = (blasint) d->n, m = (blasint) d->m;
	scs_float onef = 1.0;
	BLAS(gemv)("NoTranspose", &m, &n, &onef, A, &m, x, &one, &onef, y, &one);
}
