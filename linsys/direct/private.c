#include "private.h"

static timer linsysTimer;
static scs_float totalSolveTime;

void BLAS(potrf)(char *uplo, blasint *n, scs_float *a, blasint * lda, blasint *info);
void BLAS(trsv)(char *uplo, char *trans, char *diag, blasint *n, scs_float *a, blasint *lda, scs_float *x, blasint *incx);

char * getLinSysMethod(Data * d, Priv * p) {
	char * tmp = scs_malloc(sizeof(char) * 32);
	sprintf(tmp, "dense-direct");
	return tmp;
}

char * getLinSysSummary(Priv * p, Info * info) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "\tLin-sys: avg solve time: %1.2es\n", totalSolveTime / (info->iter + 1) / 1e3);
	totalSolveTime = 0;
	return str;
}

Priv * initPriv(Data * d) {
	blasint n = (blasint) d->n, m = (blasint) d->m, info;
	scs_float zerof = 0.0, onef = 1.0;
	scs_int k, j;
	scs_float * A = d->A->x;
	Priv * p = scs_malloc(sizeof(Priv));
	p->L = scs_calloc(n * n, sizeof(scs_float));

	/* BLAS(gemm)("Transpose", "NoTranspose", &n, &n, &m, &onef, A, &m, A, &m, &zerof, p->L, &n); */
    BLAS(syrk)("Lower", "Transpose", &n, &m, &onef, A, &m, &zerof, p->L, &n);
	for (j = 0; j < n; j++) {
		p->L[j * n + j] += d->rho_x;
	}
	BLAS(potrf)("Lower", &n, p->L, &n, &info);
	if (info != 0) {
		scs_free(p->L);
		scs_free(p);
		return NULL;
	}
	/* copy L into top half for faster solve steps */
	for (k = 0; k < n; ++k) {
		for (j = k + 1; j < n; ++j) {
			p->L[k + j * n] = p->L[j + k * n];
        }
	}
	totalSolveTime = 0.0;
	return p;
}

void freePriv(Priv * p) {
	scs_free(p->L);
	scs_free(p);
}

scs_int solveLinSys(Data * d, Priv * p, scs_float * b, const scs_float * s, scs_int iter) {
	/* returns solution to linear system */
	/* Ax = b with solution stored in b */
	scs_float * A = d->A->x;
	scs_float * L = p->L;
	blasint m = (blasint) d->m, n = (blasint) d->n, one = 1;
	scs_float onef = 1.0, negOnef = -1.0;
#ifdef EXTRAVERBOSE
	scs_printf("solving lin sys\n");
#endif
	tic(&linsysTimer);

	BLAS(gemv)("Transpose", &m, &n, &onef, A, &m, &(b[d->n]), &one, &onef, b, &one);

	/* Solve using forward-substitution, L c = b */
	BLAS(trsv)("Lower", "NoTranspose", "NonUnit", &n, L, &n, b, &one);
	/* Perform back-substitution, U x = c */
	BLAS(trsv)("Upper", "NoTranspose", "NonUnit", &n, L, &n, b, &one);

	BLAS(gemv)("NoTranspose", &m, &n, &onef, A, &m, b, &one, &negOnef, &(b[d->n]), &one);
	totalSolveTime += tocq(&linsysTimer);
#ifdef EXTRAVERBOSE
	scs_printf("finished solving lin sys\n");
#endif
	return 0;
}
