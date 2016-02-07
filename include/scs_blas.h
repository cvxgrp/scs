#ifndef SCS_BLAS_H_GUARD
#define SCS_BLAS_H_GUARD

#ifdef LAPACK_LIB_FOUND

#ifdef __cplusplus
extern "C" {
#endif

/* Default to underscore for blas / lapack */
#ifndef BLASSUFFIX
#define BLASSUFFIX _
#endif

/* annoying hack because some preprocessors can't handle empty macros */
#if defined(NOBLASSUFFIX) && NOBLASSUFFIX > 0
/* single or double precision */
#ifndef FLOAT
#define BLAS(x) d##x
#else
#define BLAS(x) s##x
#endif
#else
/* this extra indirection is needed for BLASSUFFIX to work correctly as a
 * variable */
#define stitch_(pre, x, post) pre##x##post
#define stitch__(pre, x, post) stitch_(pre, x, post)
/* single or double precision */
#ifndef FLOAT
#define BLAS(x) stitch__(d, x, BLASSUFFIX)
#else
#define BLAS(x) stitch__(s, x, BLASSUFFIX)
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

#ifdef __cplusplus
}
#endif

#endif /* LAPACK_LIB_FOUND */

#endif /* SCS_BLAS_H_GUARD */
