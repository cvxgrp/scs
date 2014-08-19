#include "private.h"

#define CG_BEST_TOL 1e-7

static idxint totCgIts;
static timer linsysTimer;
static pfloat totalSolveTime;

void BLAS(symv)(char *uplo, blasint *n, pfloat *alpha, pfloat *a, blasint *lda, const pfloat *x, blasint *incx,
		pfloat *beta, pfloat *y, blasint *incy);

pfloat BLAS(dot)(blasint *n, pfloat *dx, blasint *incx, pfloat *dy, blasint *incy);

char * getLinSysMethod(Data * d, Priv * p) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "dense-indirect, CG tol ~ 1/iter^(%2.2f)", d->cgRate);
	return str;
}

char * getLinSysSummary(Priv * p, Info * info) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
			(pfloat ) totCgIts / (info->iter + 1), totalSolveTime / (info->iter + 1) / 1e3);
	totCgIts = 0;
	totalSolveTime = 0;
	return str;
}

void getPreconditioner(Data *d, Priv *p) {
	idxint j;
	for (j = 0; j < d->n; ++j) {
		p->M[j] = 1 / p->G[j * d->n + j];
	}
}

Priv * initPriv(Data * d) {
	idxint j;
	pfloat * A = d->A->x;
	blasint n = (blasint) d->n, m = (blasint) d->m;
	pfloat onef = 1.0, zerof = 0.0;
	Priv * p = scs_malloc(sizeof(Priv));
	p->p = scs_malloc(d->n * sizeof(pfloat));
	p->z = scs_malloc(d->n * sizeof(pfloat));
	p->r = scs_malloc(d->n * sizeof(pfloat));
	p->Gp = scs_malloc(d->n * sizeof(pfloat));
	p->M = scs_malloc(d->n * sizeof(pfloat));
	/* form Gram matrix (rhoX * I+A'A) */
	p->G = scs_calloc(d->n * d->n, sizeof(pfloat));
	/* BLAS(gemm)("Transpose", "NoTranspose", &n, &n, &m, &onef, A, &m, A, &m, &zerof, p->G, &n); */
	BLAS(syrk)("Upper", "Transpose", &n, &m, &onef, A, &m, &zerof, p->G, &n);
    for (j = 0; j < d->n; j++) {
		p->G[j * d->n + j] += d->rhoX;
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

static void applyPreConditioner(pfloat * M, pfloat * z, pfloat * r, idxint n, pfloat *ipzr) {
	idxint i;
	*ipzr = 0;
	for (i = 0; i < n; ++i) {
		z[i] = r[i] * M[i];
		*ipzr += z[i] * r[i];
	}
}

/* solves (I+A'A)x = b, s warm start, solution stored in b */
static idxint pcg(Data *d, Priv * pr, const pfloat * s, pfloat * b, idxint max_its, pfloat tol) {
	idxint i = 0;
	pfloat *p = pr->p; /* cg direction */
	pfloat *Gp = pr->Gp; /* updated CG direction */
	pfloat *r = pr->r; /* cg residual */
	pfloat *G = pr->G; /* Gram matrix = (rhoX I + A'A) */
	pfloat *z = pr->z; /* preconditioned residual */
	pfloat *M = pr->M; /* preconditioner */
	pfloat ipzr, ipzrOld, alpha, negAlpha, beta;
	pfloat negOnef = -1.0, onef = 1.0, zerof = 0.0;
	blasint n = (blasint) d->n, one = 1;

	memcpy(r, b, n * sizeof(pfloat));
	if (s == NULL) {
		memset(b, 0.0, n * sizeof(pfloat));
	} else {
		BLAS(symv)("Upper", &n, &negOnef, G, &n, s, &one, &onef, r, &one);
		memcpy(b, s, n * sizeof(pfloat));
	}
	applyPreConditioner(M, z, r, n, &ipzr);
	memcpy(p, z, n * sizeof(pfloat));

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

idxint solveLinSys(Data *d, Priv * p, pfloat * b, const pfloat * s, idxint iter) {
	/* solves Mx = b, for x but stores result in b
	 s contains warm-start (if available)	p->r = b; */
	pfloat * A = d->A->x;
	idxint cgIts;
	blasint n = (blasint) d->n, m = (blasint) d->m, one = 1;
	pfloat onef = 1.0, negOnef = -1.0;
	pfloat cgTol = BLAS(nrm2)(&n, b, &one) * (iter < 0 ? CG_BEST_TOL : 1 / POWF(iter + 1, d->cgRate));
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

