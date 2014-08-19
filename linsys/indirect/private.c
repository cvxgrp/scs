#include "private.h"

#define CG_BEST_TOL 1e-7

static scs_int totCgIts;
static timer linsysTimer;
static scs_float totalSolveTime;

void BLAS(symv)(char *uplo, blasint *n, scs_float *alpha, scs_float *a, blasint *lda, const scs_float *x, blasint *incx,
		scs_float *beta, scs_float *y, blasint *incy);

scs_float BLAS(dot)(blasint *n, scs_float *dx, blasint *incx, scs_float *dy, blasint *incy);

char * getLinSysMethod(Data * d, Priv * p) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "dense-indirect, CG tol ~ 1/iter^(%2.2f)", d->cg_rate);
	return str;
}

char * getLinSysSummary(Priv * p, Info * info) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
			(scs_float ) totCgIts / (info->iter + 1), totalSolveTime / (info->iter + 1) / 1e3);
	totCgIts = 0;
	totalSolveTime = 0;
	return str;
}

void getPreconditioner(Data *d, Priv *p) {
	scs_int j;
	for (j = 0; j < d->n; ++j) {
		p->M[j] = 1 / p->G[j * d->n + j];
	}
}

Priv * initPriv(Data * d) {
	scs_int j;
	scs_float * A = d->A->x;
	blasint n = (blasint) d->n, m = (blasint) d->m;
	scs_float onef = 1.0, zerof = 0.0;
	Priv * p = scs_malloc(sizeof(Priv));
	p->p = scs_malloc(d->n * sizeof(scs_float));
	p->z = scs_malloc(d->n * sizeof(scs_float));
	p->r = scs_malloc(d->n * sizeof(scs_float));
	p->Gp = scs_malloc(d->n * sizeof(scs_float));
	p->M = scs_malloc(d->n * sizeof(scs_float));
	/* form Gram matrix (rho_x * I+A'A) */
	p->G = scs_calloc(d->n * d->n, sizeof(scs_float));
	/* BLAS(gemm)("Transpose", "NoTranspose", &n, &n, &m, &onef, A, &m, A, &m, &zerof, p->G, &n); */
	BLAS(syrk)("Upper", "Transpose", &n, &m, &onef, A, &m, &zerof, p->G, &n);
    for (j = 0; j < d->n; j++) {
		p->G[j * d->n + j] += d->rho_x;
	}
	getPreconditioner(d, p);
	totalSolveTime = 0;
	totCgIts = 0;
	return p;
}

void freePriv(Priv * p) {
	if (p) {
		if (p->p)
			scs_free(p->p);
		if (p->r)
			scs_free(p->r);
		if (p->Gp)
			scs_free(p->Gp);
		if (p->z)
			scs_free(p->z);
		if (p->M)
			scs_free(p->M);
		if (p->G)
			scs_free(p->G);
		scs_free(p);
	}
}

static void applyPreConditioner(scs_float * M, scs_float * z, scs_float * r, scs_int n, scs_float *ipzr) {
	scs_int i;
	*ipzr = 0;
	for (i = 0; i < n; ++i) {
		z[i] = r[i] * M[i];
		*ipzr += z[i] * r[i];
	}
}

/* solves (I+A'A)x = b, s warm start, solution stored in b */
static scs_int pcg(Data *d, Priv * pr, const scs_float * s, scs_float * b, scs_int max_its, scs_float tol) {
	scs_int i = 0;
	scs_float *p = pr->p; /* cg direction */
	scs_float *Gp = pr->Gp; /* updated CG direction */
	scs_float *r = pr->r; /* cg residual */
	scs_float *G = pr->G; /* Gram matrix = (rho_x I + A'A) */
	scs_float *z = pr->z; /* preconditioned residual */
	scs_float *M = pr->M; /* preconditioner */
	scs_float ipzr, ipzrOld, alpha, negAlpha, beta;
	scs_float negOnef = -1.0, onef = 1.0, zerof = 0.0;
	blasint n = (blasint) d->n, one = 1;

	memcpy(r, b, n * sizeof(scs_float));
	if (s == NULL) {
		memset(b, 0.0, n * sizeof(scs_float));
	} else {
		BLAS(symv)("Upper", &n, &negOnef, G, &n, s, &one, &onef, r, &one);
		memcpy(b, s, n * sizeof(scs_float));
	}
	applyPreConditioner(M, z, r, n, &ipzr);
	memcpy(p, z, n * sizeof(scs_float));

	for (i = 0; i < max_its; i++) {
		BLAS(symv)("Upper", &n, &onef, G, &n, p, &one, &zerof, Gp, &one);

		alpha = ipzr / BLAS(dot)(&n, p, &one, Gp, &one);
		negAlpha = -alpha;

		BLAS(axpy)(&n, &alpha, p, &one, b, &one);
		BLAS(axpy)(&n, &negAlpha, Gp, &one, r, &one);

		if (BLAS(nrm2)(&n, r, &one) < tol) {
			return i + 1;
		}
		ipzrOld = ipzr;
		applyPreConditioner(M, z, r, n, &ipzr);

		beta = ipzr / ipzrOld;
		BLAS(scal)(&n, &beta, p, &one);
		BLAS(axpy)(&n, &onef, z, &one, p, &one);
	}
	return i;
}

scs_int solveLinSys(Data *d, Priv * p, scs_float * b, const scs_float * s, scs_int iter) {
	/* solves Mx = b, for x but stores result in b
	 s contains warm-start (if available)	p->r = b; */
	scs_float * A = d->A->x;
	scs_int cgIts;
	blasint n = (blasint) d->n, m = (blasint) d->m, one = 1;
	scs_float onef = 1.0, negOnef = -1.0;
	scs_float cgTol = BLAS(nrm2)(&n, b, &one) * (iter < 0 ? CG_BEST_TOL : 1 / POWF(iter + 1, d->cg_rate));
#ifdef EXTRAVERBOSE
	scs_printf("solving lin sys\n");
#endif
	tic(&linsysTimer);
	BLAS(gemv)("Transpose", &m, &n, &onef, A, &m, &(b[d->n]), &one, &onef, b, &one);
	cgIts = pcg(d, p, s, b, d->n, MAX(cgTol, CG_BEST_TOL));
	BLAS(gemv)("NoTranpose", &m, &n, &onef, A, &m, b, &one, &negOnef, &(b[d->n]), &one);
#ifdef EXTRAVERBOSE
	scs_printf("\tCG iterations: %i\n", (int) cgIts);
#endif
	if (iter >= 0) {
		totCgIts += cgIts;
	}

	totalSolveTime += tocq(&linsysTimer);
	return 0;
}

