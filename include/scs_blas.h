#ifndef SCS_BLAS_H_GUARD
#define SCS_BLAS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

/* check for empty, indirection required for stiching macros */
#define _SCS_EMPTY_HELPER
#define EMPTY(name) defined(name ## _SCS_EMPTY_HELPER)
#define IS_EMPTY(name) EMPTY(name)

/* Default to underscore for blas / lapack */
#ifndef BLASSUFFIX
#define BLASSUFFIX _
#endif

#ifdef LAPACK_LIB_FOUND
    /* empty BLASSUFFIX macros are causing errors and warnings */
    #if !IS_EMPTY(BLASSUFFIX)
        /* this extra indirection is needed for BLASSUFFIX to work correctly as a variable */
        #define stitch_(pre,x,post) pre ## x ## post
        #define stitch__(pre,x,post) stitch_(pre,x,post)
        /* single or double precision */
        #ifndef FLOAT
            #define BLAS(x) stitch__(d,x,BLASSUFFIX)
        #else
            #define BLAS(x) stitch__(s,x,BLASSUFFIX)
        #endif
    #else
        #ifndef FLOAT
            #define BLAS(x) a ## x
        #else
            #define BLAS(x) s ## x
        #endif
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
#endif
