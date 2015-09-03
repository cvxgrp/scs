#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include "amatrix.h"

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

/* Default to underscore for blas / lapack */
#ifndef BLASSUFFIX
#define BLASSUFFIX _
#endif

/* this extra indirection is needed for BLASSUFFIX to work correctly as a variable */
#define stitch_(pre,x,post) pre ## x ## post
#define stitch__(pre,x,post) stitch_(pre,x,post)
/* single or double precision */
#ifndef FLOAT
#define BLAS(x) stitch__(d,x,BLASSUFFIX)
#else
#define BLAS(x) stitch__(s,x,BLASSUFFIX)
#endif

#ifdef MATLAB_MEX_FILE
typedef ptrdiff_t blasint;
#elif defined BLAS64
#include <stdint.h>
typedef int64_t blasint;
#else
typedef int blasint;
#endif

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
