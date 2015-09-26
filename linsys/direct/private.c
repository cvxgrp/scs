#include "private.h"

static timer linsysTimer;
static scs_float totalSolveTime;

void BLAS(potrf)(char *uplo, blasint *n, scs_float *a, blasint * lda, blasint *info);
void BLAS(trsv)(char *uplo, char *trans, char *diag, blasint *n, scs_float *a, blasint *lda, scs_float *x, blasint *incx);

char * getLinSysMethod(const AMatrix * A, const Settings * stgs) {
	char * tmp = scs_malloc(sizeof(char) * 32);
	sprintf(tmp, "dense-direct");
	return tmp;
}

char * getLinSysSummary(Priv * p, const Info * info) {
	char * str = scs_malloc(sizeof(char) * 128);
	sprintf(str, "\tLin-sys: avg solve time: %1.2es\n", totalSolveTime / (info->iter + 1) / 1e3);
	totalSolveTime = 0;
	return str;
}

Priv * initPriv(const AMatrix * A, const Settings * stgs) {
	blasint n = (blasint) A->n, m = (blasint) A->m, info;
	scs_float zerof = 0.0, onef = 1.0;
	scs_int k, j;
	Priv * p = scs_malloc(sizeof(Priv));
	p->L = scs_calloc(n * n, sizeof(scs_float));

    BLAS(syrk)("Lower", "Transpose", &n, &m, &onef, A->x, &m, &zerof, p->L, &n);
	for (j = 0; j < n; j++) {
		p->L[j * n + j] += stgs->rho_x;
	}
	BLAS(potrf)("Lower", &n, p->L, &n, &info);
	if (info != 0) {
		scs_free(p->L);
		scs_free(p);
		return SCS_NULL;
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

scs_int solveLinSys(const AMatrix * A, const Settings * stgs, Priv * p, scs_float * b, const scs_float * s, scs_int iter) {
	/* returns solution to linear system */
	/* Ax = b with solution stored in b */
	scs_float * L = p->L;
	blasint m = (blasint) A->m, n = (blasint) A->n, one = 1;
	scs_float onef = 1.0, negOnef = -1.0;
#ifdef EXTRAVERBOSE
	scs_printf("solving lin sys\n");
#endif
	tic(&linsysTimer);

	BLAS(gemv)("Transpose", &m, &n, &onef, A->x, &m, &(b[A->n]), &one, &onef, b, &one);

	/* Solve using forward-substitution, L c = b */
	BLAS(trsv)("Lower", "NoTranspose", "NonUnit", &n, L, &n, b, &one);
	/* Perform back-substitution, U x = c */
	BLAS(trsv)("Upper", "NoTranspose", "NonUnit", &n, L, &n, b, &one);

	BLAS(gemv)("NoTranspose", &m, &n, &onef, A->x, &m, b, &one, &negOnef, &(b[A->n]), &one);
	totalSolveTime += tocq(&linsysTimer);
#ifdef EXTRAVERBOSE
	scs_printf("finished solving lin sys\n");
#endif
	return 0;
}

