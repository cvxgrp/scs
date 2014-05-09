#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#include "cones.h"
#include "amatrix.h"

/* contains routines common to direct and indirect sparse solvers */

#define MIN_SCALE 1e-2
#define MAX_SCALE 1e3

idxint validateLinSys(Data *d) {
	AMatrix * A = d->A;
	idxint i, rMax, Anz;
	if (!A->x || !A->i || !A->p) {
		scs_printf("data incompletely specified\n");
		return -1;
	}
	/* detects degenerate problems, sometimes not wanted:
	 for (i = 0; i < d->n; ++i) {
	 	 if (A->p[i] >= A->p[i + 1]) {
	 	 	 scs_printf("A->p not strictly increasing\n");
	 	 	 return -1;
	 	 }
	 }
	*/
	Anz = A->p[d->n];
	if (((pfloat) Anz / d->m > d->n) || (Anz <= 0)) {
		scs_printf("Anz (nonzeros in A) = %i, outside of valid range\n", (int) Anz);
		return -1;
	}
	rMax = 0;
	for (i = 0; i < Anz; ++i) {
		if (A->i[i] > rMax)
			rMax = A->i[i];
	}
	if (rMax > d->m - 1) {
		scs_printf("number of rows in A inconsistent with input dimension\n");
		return -1;
	}
	return 0;
}

void normalizeA(Data * d, Work * w, Cone * k) {
	AMatrix * A = d->A;
	pfloat * D = scs_calloc(d->m, sizeof(pfloat));
	pfloat * E = scs_calloc(d->n, sizeof(pfloat));
	idxint i, j, count, delta, *boundaries, c1, c2;
	pfloat wrk, *nms, e;
	idxint numBoundaries = getConeBoundaries(k, &boundaries);
#ifdef EXTRAVERBOSE
	timer normalizeTimer;
	scs_printf("normalizing A\n");
	tic(&normalizeTimer);
#endif

	/* calculate row norms */
	for (i = 0; i < d->n; ++i) {
		c1 = A->p[i];
		c2 = A->p[i + 1];
		for (j = c1; j < c2; ++j) {
			wrk = A->x[j];
			D[A->i[j]] += wrk * wrk;
		}
	}
	for (i = 0; i < d->m; ++i) {
		D[i] = sqrt(D[i]); /* just the norms */
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
		if (D[i] < MIN_SCALE)
			D[i] = 1;
		else if (D[i] > MAX_SCALE)
			D[i] = MAX_SCALE;
	}

	/* scale the rows with D */
	for (i = 0; i < d->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			A->x[j] /= D[A->i[j]];
		}
	}
	/* calculate and scale by col norms, E */
	for (i = 0; i < d->n; ++i) {
		c1 = A->p[i + 1] - A->p[i];
		e = calcNorm(&(A->x[A->p[i]]), c1);
		if (e < MIN_SCALE)
			e = 1;
		else if (e > MAX_SCALE)
			e = MAX_SCALE;
		scaleArray(&(A->x[A->p[i]]), 1.0 / e, c1);
		E[i] = e;
	}

	nms = scs_calloc(d->m, sizeof(pfloat));
	for (i = 0; i < d->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			wrk = A->x[j];
			nms[A->i[j]] += wrk * wrk;
		}
	}
	w->meanNormRowA = 0.0;
	for (i = 0; i < d->m; ++i) {
		w->meanNormRowA += sqrt(nms[i]) / d->m;
	}
	scs_free(nms);

	if (d->SCALE != 1) {
		scaleArray(A->x, d->SCALE, A->p[d->n]);
	}

	w->D = D;
	w->E = E;

#ifdef EXTRAVERBOSE
	scs_printf("finished normalizing A, time: %6f s\n", tocq(&normalizeTimer) / 1e3);
#endif
}

void unNormalizeA(Data *d, Work * w) {
	idxint i, j;
	pfloat * D = w->D;
	pfloat * E = w->E;
	AMatrix * A = d->A;
	for (i = 0; i < d->n; ++i) {
		scaleArray(&(A->x[A->p[i]]), E[i] / d->SCALE, A->p[i + 1] - A->p[i]);
	}
	for (i = 0; i < d->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			A->x[j] *= D[A->i[j]];
		}
	}
}

#endif
