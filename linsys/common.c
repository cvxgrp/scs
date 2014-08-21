#include "common.h"
/* contains routines common to direct and indirect sparse solvers */

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)
#define NUM_SCALE_PASSES 1 /* additional passes don't help much */

scs_int validateLinSys(Data *d) {
	AMatrix * A = d->A;
	scs_int i, rMax, Anz;
	if (!A->x || !A->i || !A->p) {
		scs_printf("data incompletely specified\n");
		return -1;
	}
	/* detects degenerate problems, typically this check is not wanted:
	 for (i = 0; i < d->n; ++i) {
	 if (A->p[i] >= A->p[i + 1]) {
	 scs_printf("A->p not strictly increasing\n");
	 return -1;
	 }
	 }
	 */
	Anz = A->p[d->n];
	if (((scs_float) Anz / d->m > d->n) || (Anz <= 0)) {
		scs_printf("Anz (nonzeros in A) = %li, outside of valid range\n", (long) Anz);
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

void printAMatrix(Data * d) {
	scs_int i, j;
	AMatrix * A = d->A;
	/* TODO: this is to prevent clogging stdout */
	if (A->p[d->n] < 2500) {
		scs_printf("\n");
		for (i = 0; i < d->n; ++i) {
			scs_printf("Col %li: ", (long) i);
			for (j = A->p[i]; j < A->p[i + 1]; j++) {
				scs_printf("A[%li,%li] = %4f, ", (long) A->i[j], (long) i, A->x[j]);
			}
			scs_printf("norm col = %4f\n", calcNorm(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
		}
		scs_printf("norm A = %4f\n", calcNorm(A->x, A->p[d->n]));
	}
}

void normalizeA(Data * d, Work * w, Cone * k) {
	AMatrix * A = d->A;
	scs_float * D = scs_malloc(d->m * sizeof(scs_float));
	scs_float * E = scs_malloc(d->n * sizeof(scs_float));
	scs_float * Dt = scs_malloc(d->m * sizeof(scs_float));
	scs_float * Et = scs_malloc(d->n * sizeof(scs_float));
	scs_float * nms = scs_calloc(d->m, sizeof(scs_float));
	scs_float minRowScale = MIN_SCALE * SQRTF((scs_float) d->n), maxRowScale = MAX_SCALE * SQRTF((scs_float) d->n);
	scs_float minColScale = MIN_SCALE * SQRTF((scs_float) d->m), maxColScale = MAX_SCALE * SQRTF((scs_float) d->m);
	scs_int i, j, l, count, delta, *boundaries, c1, c2;
	scs_float wrk, e;
	scs_int numBoundaries = getConeBoundaries(k, &boundaries);

#ifdef EXTRAVERBOSE
	timer normalizeTimer;
	tic(&normalizeTimer);
	scs_printf("normalizing A\n");
	printAMatrix(d);
#endif

	for (l = 0; l < NUM_SCALE_PASSES; ++l) {
		memset(D, 0, d->m * sizeof(scs_float));
		memset(E, 0, d->n * sizeof(scs_float));
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

		for (i = 0; i < d->m; ++i) {
			if (D[i] < minRowScale)
				D[i] = 1;
			else if (D[i] > maxRowScale)
				D[i] = maxRowScale;
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
			if (e < minColScale)
				e = 1;
			else if (e > maxColScale)
				e = maxColScale;
			scaleArray(&(A->x[A->p[i]]), 1.0 / e, c1);
			E[i] = e;
		}

		for (i = 0; i < d->m; ++i) {
			Dt[i] = (l == 0) ? D[i] : Dt[i] * D[i];
		}
		for (i = 0; i < d->n; ++i) {
			Et[i] = (l == 0) ? E[i] : Et[i] * E[i];
		}
	}
	scs_free(boundaries);
	scs_free(D);
	scs_free(E);

	/* calculate mean of row norms of A */
	for (i = 0; i < d->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			wrk = A->x[j];
			nms[A->i[j]] += wrk * wrk;
		}
	}
	w->meanNormRowA = 0.0;
	for (i = 0; i < d->m; ++i) {
		w->meanNormRowA += SQRTF(nms[i]) / d->m;
	}
	scs_free(nms);

	/* calculate mean of col norms of A */
	w->meanNormColA = 0.0;
	for (i = 0; i < d->n; ++i) {
		c1 = A->p[i + 1] - A->p[i];
		w->meanNormColA += calcNorm(&(A->x[A->p[i]]), c1) / d->n;
	}

	/* scale up by d->SCALE if not equal to 1 */
	if (d->scale != 1) {
		scaleArray(A->x, d->scale, A->p[d->n]);
	}

	w->D = Dt;
	w->E = Et;

#ifdef EXTRAVERBOSE
	scs_printf("finished normalizing A, time: %1.2es\n", tocq(&normalizeTimer) / 1e3);
	printAMatrix(d);
#endif
}

void unNormalizeA(Data *d, Work * w) {
	scs_int i, j;
	scs_float * D = w->D;
	scs_float * E = w->E;
	AMatrix * A = d->A;
	for (i = 0; i < d->n; ++i) {
		scaleArray(&(A->x[A->p[i]]), E[i] / d->scale, A->p[i + 1] - A->p[i]);
	}
	for (i = 0; i < d->n; ++i) {
		for (j = A->p[i]; j < A->p[i + 1]; ++j) {
			A->x[j] *= D[A->i[j]];
		}
	}
}
