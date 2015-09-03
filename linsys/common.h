#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include "amatrix.h"
#include "scs_blas.h"

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

void BLAS(scal)(blasint *n, const scs_float *da, scs_float *dx, blasint *incx);
scs_float BLAS(nrm2)(blasint *n, scs_float *x, blasint *incx);
void BLAS(gemv)(char *trans, blasint *m, blasint *n, scs_float * alpha, scs_float *a, blasint *lda, const scs_float *x, blasint *incx,
		scs_float *beta, scs_float *y, blasint *incy);
void BLAS(axpy)(blasint *n, scs_float *da, scs_float *dx, blasint *incx, scs_float *dy, blasint *incy_);
void BLAS(syrk)(char *uplo, char *trans, blasint *n, blasint *k, scs_float *alpha, scs_float *a, blasint *lda, scs_float *beta, 
        scs_float *c, blasint *ldc);

/* void BLAS(gemm)(char *transa, char *transb, blasint *m, blasint * n, blasint *k, scs_float *alpha, scs_float *a, blasint *lda,
        scs_float *b, blasint *ldb, scs_float *beta, scs_float *c, blasint *ldc); */

#ifdef __cplusplus
}
#endif

#endif
