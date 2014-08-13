#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#include "glbopts.h"
#include "cones.h"
#include "amatrix.h"

#define MIN_SCALE 1e-3
#define MAX_SCALE 1e3

/* underscore for blas / lapack, single or double precision */
#ifdef NOBLASUNDERSCORE
	#ifndef FLOAT
		#define BLAS(x) d ## x
	#else
		#define BLAS(x) s ## x
	#endif
#else
	#ifndef FLOAT
		#define BLAS(x) d ## x ## _
	#else
		#define BLAS(x) s ## x ## _
	#endif
#endif

#ifdef MATLAB_MEX_FILE
typedef ptrdiff_t blasint;
#elif defined BLAS64
#include <stdint.h>
typedef int64_t blasint;
#else
typedef int blasint;
#endif

void BLAS(scal)(blasint *n, pfloat *da, pfloat *dx, blasint *incx);
pfloat BLAS(nrm2)(blasint *n, pfloat *x, blasint *incx);
void BLAS(gemv)(char *trans, blasint *m, blasint *n, pfloat * alpha, pfloat *a, blasint *lda, const pfloat *x, blasint *incx,
		pfloat *beta, pfloat *y, blasint *incy);
void BLAS(axpy)(blasint *n, pfloat *da, pfloat *dx, blasint *incx, pfloat *dy, blasint *incy_);
void BLAS(syrk)(char *uplo, char *trans, blasint *n, blasint *k, pfloat *alpha, pfloat *a, blasint *lda, pfloat *beta, 
        pfloat *c, blasint *ldc);

/* void BLAS(gemm)(char *transa, char *transb, blasint *m, blasint * n, blasint *k, pfloat *alpha, pfloat *a, blasint *lda,
        pfloat *b, blasint *ldb, pfloat *beta, pfloat *c, blasint *ldc); */

idxint validateLinSys(Data *d);
void normalizeA(Data * d, Work * w, Cone * k);
void unNormalizeA(Data *d, Work * w);
void accumByAtrans(Data * d, Priv * p, const pfloat * x, pfloat * y);
#endif
