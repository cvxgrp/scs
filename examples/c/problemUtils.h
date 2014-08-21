#ifndef PUTILS_H_GUARD
#define PUTILS_H_GUARD

#include "scs.h"
#include "linsys/amatrix.h"
#include "linsys/common.h"

#define PI (3.141592654)
#ifdef DLONG
#ifdef _WIN64
/* this is a Microsoft extension, but also works with MinGW-w64 */
#define INTRW "%I64d"
#else
#define INTRW "%ld"
#endif
#else
#define INTRW "%i"
#endif

#ifndef FLOAT
#define FLOATRW "%lf"
#else
#define FLOATRW "%f"
#endif

/* uniform random number in [-1,1] */
scs_float rand_scs_float(void) {
	return (2 * (((scs_float) rand()) / RAND_MAX) - 1);
}

/* normal random var */
static scs_float U, V;
static scs_int phase = 0;
scs_float rand_gauss(void) {
	scs_float Z;
	if (phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = SQRTF(-2 * log(U)) * sin(2 * PI * V);
	} else
		Z = SQRTF(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;
	return Z;
}

void perturbVector(scs_float * v, scs_int l) {
	scs_int i;
	for (i = 0; i < l; i++) {
		v[i] += 0.01 * rand_gauss();
	}
}

void freeData(Data * d, Cone * k) {
	if (d) {
		if (d->b)
			scs_free(d->b);
		if (d->c)
			scs_free(d->c);
		if (d->A) {
			if (d->A->x)
				scs_free(d->A->x);
			scs_free(d->A);
		}
		scs_free(d);
	}
	if (k) {
		if (k->q)
			scs_free(k->q);
		if (k->s)
			scs_free(k->s);
		scs_free(k);
	}
	d = NULL;
	k = NULL;
}

void freeSol(Sol *sol) {
	if (sol) {
		if (sol->x) {
			scs_free(sol->x);
			sol->x = NULL;
		}
		if (sol->y) {
			scs_free(sol->y);
			sol->y = NULL;
		}
		if (sol->s) {
			scs_free(sol->s);
			sol->s = NULL;
		}
        scs_free(sol);
	}
    sol = NULL;
}

void genRandomProbData(Data * d, Cone * k, Sol * opt_sol) {
	blasint one = 1, n = (blasint) d->n, m = (blasint) d->m;
	scs_float onef = 1.0, negOnef = -1.0, zerof = 0.0;
	AMatrix * A = d->A = scs_calloc(1, sizeof(AMatrix));
	scs_float * b = d->b = scs_calloc(m, sizeof(scs_float));
	scs_float * c = d->c = scs_calloc(n, sizeof(scs_float));
	scs_float * x = opt_sol->x = scs_calloc(n, sizeof(scs_float));
	scs_float * y = opt_sol->y = scs_calloc(m, sizeof(scs_float));
	scs_float * s = opt_sol->s = scs_calloc(m, sizeof(scs_float));
	/* temporary variables */
	scs_float * z = scs_calloc(m, sizeof(scs_float));
	scs_int i;

	A->x = scs_calloc(n * m, sizeof(scs_float));
	/* y, s >= 0 and y'*s = 0 */
	for (i = 0; i < m; i++) {
		y[i] = z[i] = rand_scs_float();
	}

	projDualCone(y, k, NULL, -1);

	for (i = 0; i < m; i++) {
		b[i] = s[i] = y[i] - z[i];
	}

	for (i = 0; i < n; i++) {
		x[i] = rand_scs_float();
	}

	/*
	 * c = -A'*y
	 * b = A*x + s
	 */
    for (i = 0; i < n * m; i++) {
        A->x[i] = rand_scs_float();
    }
	BLAS(gemv)("NoTranspose", &m, &n, &onef, A->x, &m, x, &one, &onef, b, &one);
	BLAS(gemv)("Transpose", &m, &n, &negOnef, A->x, &m, y, &one, &zerof, c, &one);
    scs_printf("Finished generating random cone prob\n");
	scs_free(z);
}

#endif
