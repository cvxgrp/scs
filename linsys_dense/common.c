#include "common.h"

idxint validateLinSys(Data *d) {
	return (d->A && d->A->x);
}

void normalizeA(Data * d, Work * w, Cone * k) {
	idxint i, j, count, delta;
	pfloat wrk;
    pfloat * D = scs_calloc(d->m, sizeof(pfloat));
    pfloat * E = scs_calloc(d->n, sizeof(pfloat));
    pfloat minRowScale = MIN_SCALE * SQRTF(d->n), maxRowScale = MAX_SCALE * SQRTF(d->n);
    pfloat minColScale = MIN_SCALE * SQRTF(d->m), maxColScale = MAX_SCALE * SQRTF(d->m);
    idxint *boundaries, numBoundaries = getConeBoundaries(k, &boundaries);
    blasint n = (blasint) d->n, m = (blasint) d->m, one = 1, nm = n * m;

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

	if (d->SCALE != 1) {
		BLAS(scal)(&nm, &d->SCALE, d->A->x, &one);
	}

	w->D = D;
	w->E = E;
}

void unNormalizeA(Data *d, Work * w) {
	idxint i;
	pfloat * D = w->D;
	pfloat * E = w->E;
	pfloat invScale = 1.0 / d->SCALE;
	blasint n = (blasint) d->n, m = (blasint) d->m, one = 1, nm = n * m;
	for (i = 0; i < d->n; ++i) {
		BLAS(scal)(&m, &(E[i]), &(d->A->x[i * d->m]), &one);
	}
	for (i = 0; i < d->m; ++i) {
		BLAS(scal)(&n, &(D[i]), &(d->A->x[i]), &m);
	}
	BLAS(scal)(&nm, &invScale, d->A->x, &one);
}

void accumByAtrans(Data * d, Priv * p, const pfloat * x, pfloat * y) {
	pfloat * A = d->A->x;
	blasint one = 1, n = (blasint) d->n, m = (blasint) d->m;
	pfloat onef = 1.0;
	BLAS(gemv)("Transpose", &m, &n, &onef, A, &m, x, &one, &onef, y, &one);
}
void accumByA(Data * d, Priv * p, const pfloat * x, pfloat * y) {
	pfloat * A = d->A->x;
	blasint one = 1, n = (blasint) d->n, m = (blasint) d->m;
	pfloat onef = 1.0;
	BLAS(gemv)("NoTranspose", &m, &n, &onef, A, &m, x, &one, &onef, y, &one);
}
