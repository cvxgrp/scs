#include "private.h"

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1
#define PRINT_INTERVAL 100

static scs_int totCgIts;
static timer linsysTimer;
static scs_float totalSolveTime;

char * getLinSysMethod(Data * d, Priv * p) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "sparse-indirect, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)", (long ) d->A->p[d->n], d->cg_rate);
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

/* M = inv ( diag ( RHO_X * I + A'A ) ) */
void getPreconditioner(Data *d, Priv *p) {
	scs_int i;
	scs_float * M = p->M;
	AMatrix * A = d->A;

#ifdef EXTRAVERBOSE
	scs_printf("getting pre-conditioner\n");
#endif

	for (i = 0; i < d->n; ++i) {
		M[i] = 1 / (d->rho_x + calcNormSq(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
		/* M[i] = 1; */
	}

#ifdef EXTRAVERBOSE
	scs_printf("finished getting pre-conditioner\n");
#endif

}

static void transpose(Data * d, Priv * p) {
	scs_int * Ci = p->Ati;
	scs_int * Cp = p->Atp;
	scs_float * Cx = p->Atx;
	scs_int m = d->m;
	scs_int n = d->n;

	scs_int * Ap = d->A->p;
	scs_int * Ai = d->A->i;
	scs_float * Ax = d->A->x;

	scs_int i, j, q, *z, c1, c2;
#ifdef EXTRAVERBOSE
	timer transposeTimer;
	scs_printf("transposing A\n");
	tic(&transposeTimer);
#endif

	z = scs_calloc(m, sizeof(scs_int));
	for (i = 0; i < Ap[n]; i++)
		z[Ai[i]]++; /* row counts */
	cs_cumsum(Cp, z, m); /* row pointers */

	for (j = 0; j < n; j++) {
		c1 = Ap[j];
		c2 = Ap[j + 1];
		for (i = c1; i < c2; i++) {
			q = z[Ai[i]];
			Ci[q] = j; /* place A(i,j) as entry C(j,i) */
			Cx[q] = Ax[i];
			z[Ai[i]]++;
		}
	}
	scs_free(z);

#ifdef EXTRAVERBOSE
	scs_printf("finished transposing A, time: %1.2es\n", tocq(&transposeTimer) / 1e3);
#endif

}

void freePriv(Priv * p) {
	if (p) {
		if (p->p)
			scs_free(p->p);
		if (p->r)
			scs_free(p->r);
		if (p->Gp)
			scs_free(p->Gp);
		if (p->tmp)
			scs_free(p->tmp);
		if (p->Ati)
			scs_free(p->Ati);
		if (p->Atx)
			scs_free(p->Atx);
		if (p->Atp)
			scs_free(p->Atp);
		if (p->z)
			scs_free(p->z);
		if (p->M)
			scs_free(p->M);
		scs_free(p);
	}
}

/* solves (I+A'A)x = b, s warm start, solution stored in b */
/*y = (RHO_X * I + A'A)x */
static void matVec(Data * d, Priv * p, const scs_float * x, scs_float * y) {
	scs_float * tmp = p->tmp;
	memset(tmp, 0, d->m * sizeof(scs_float));
	accumByA(d, p, x, tmp);
	memset(y, 0, d->n * sizeof(scs_float));
	accumByAtrans(d, p, tmp, y);
	addScaledArray(y, x, d->n, d->rho_x);
}

void _accumByAtrans(scs_int n, scs_float * Ax, scs_int * Ai, scs_int * Ap, const scs_float *x, scs_float *y) {
	/* y  = A'*x
	 A in column compressed format
	 parallelizes over columns (rows of A')
	 */
	scs_int p, j;
	scs_int c1, c2;
	scs_float yj;
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
}

void accumByAtrans(Data * d, Priv * p, const scs_float *x, scs_float *y) {
	AMatrix * A = d->A;
	_accumByAtrans(d->n, A->x, A->i, A->p, x, y);
}
void accumByA(Data * d, Priv * p, const scs_float *x, scs_float *y) {
	_accumByAtrans(d->m, p->Atx, p->Ati, p->Atp, x, y);
}
static void applyPreConditioner(scs_float * M, scs_float * z, scs_float * r, scs_int n, scs_float *ipzr) {
	scs_int i;
	*ipzr = 0;
	for (i = 0; i < n; ++i) {
		z[i] = r[i] * M[i];
		*ipzr += z[i] * r[i];
	}
}

Priv * initPriv(Data * d) {
	AMatrix * A = d->A;
	Priv * p = scs_calloc(1, sizeof(Priv));
	p->p = scs_malloc((d->n) * sizeof(scs_float));
	p->r = scs_malloc((d->n) * sizeof(scs_float));
	p->Gp = scs_malloc((d->n) * sizeof(scs_float));
	p->tmp = scs_malloc((d->m) * sizeof(scs_float));

	/* preconditioner memory */
	p->z = scs_malloc((d->n) * sizeof(scs_float));
	p->M = scs_malloc((d->n) * sizeof(scs_float));

	p->Ati = scs_malloc((A->p[d->n]) * sizeof(scs_int));
	p->Atp = scs_malloc((d->m + 1) * sizeof(scs_int));
	p->Atx = scs_malloc((A->p[d->n]) * sizeof(scs_float));
	transpose(d, p);
	getPreconditioner(d, p);
	totalSolveTime = 0;
	totCgIts = 0;
	if (!p->p || !p->r || !p->Gp || !p->tmp || !p->Ati || !p->Atp || !p->Atx) {
		freePriv(p);
		return NULL;
	}
	return p;
}

static scs_int pcg(Data *d, Priv * pr, const scs_float * s, scs_float * b, scs_int max_its, scs_float tol) {
	scs_int i, n = d->n;
	scs_float ipzr, ipzrOld, alpha;
	scs_float *p = pr->p; /* cg direction */
	scs_float *Gp = pr->Gp; /* updated CG direction */
	scs_float *r = pr->r; /* cg residual */
	scs_float *z = pr->z; /* for preconditioning */
	scs_float *M = pr->M; /* inverse diagonal preconditioner */

	if (s == NULL) {
		memcpy(r, b, n * sizeof(scs_float));
		memset(b, 0, n * sizeof(scs_float));
	} else {
		matVec(d, pr, s, r);
		addScaledArray(r, b, n, -1);
		scaleArray(r, -1, n);
		memcpy(b, s, n * sizeof(scs_float));
	}
	applyPreConditioner(M, z, r, n, &ipzr);
	memcpy(p, z, n * sizeof(scs_float));

	for (i = 0; i < max_its; ++i) {
		matVec(d, pr, p, Gp);

		alpha = ipzr / innerProd(p, Gp, n);
		addScaledArray(b, p, n, alpha);
		addScaledArray(r, Gp, n, -alpha);

		if (calcNorm(r, n) < tol) {
            #ifdef EXTRAVERBOSE
            scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, calcNorm(r, n), (long) i+1);
            #endif
			return i + 1;
		}
		ipzrOld = ipzr;
		applyPreConditioner(M, z, r, n, &ipzr);

		scaleArray(p, ipzr / ipzrOld, n);
		addScaledArray(p, z, n, 1);
	}
	return i;
}

scs_int solveLinSys(Data *d, Priv * p, scs_float * b, const scs_float * s, scs_int iter) {
	scs_int cgIts;
	scs_float cgTol = calcNorm(b, d->n) * (iter < 0 ? CG_BEST_TOL : CG_MIN_TOL / POWF((scs_float) iter + 1, d->cg_rate));

	tic(&linsysTimer);
	/* solves Mx = b, for x but stores result in b */
	/* s contains warm-start (if available) */
	accumByAtrans(d, p, &(b[d->n]), b);
	/* solves (I+A'A)x = b, s warm start, solution stored in b */
	cgIts = pcg(d, p, s, b, d->n, MAX(cgTol, CG_BEST_TOL));
	scaleArray(&(b[d->n]), -1, d->m);
	accumByA(d, p, b, &(b[d->n]));

	if (iter >= 0) {
		totCgIts += cgIts;
	}

	totalSolveTime += tocq(&linsysTimer);
#ifdef EXTRAVERBOSE
	scs_printf("linsys solve time: %1.2es\n", tocq(&linsysTimer) / 1e3);
#endif
	return 0;
}

