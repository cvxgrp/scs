#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#include <math.h>

/* redefine printfs and memory allocators as needed */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define scs_printf   mexPrintf
#define scs_free     mxFree
#define scs_malloc   mxMalloc
#define scs_calloc   mxCalloc
#elif defined PYTHON
#include <Python.h>
#include <stdlib.h>
#define scs_printf   PySys_WriteStdout
#define scs_free     free
#define scs_malloc   malloc
#define scs_calloc   calloc
#else
#include <stdio.h>
#include <stdlib.h>
#define scs_printf   printf
#define scs_free     free
#define scs_malloc   malloc
#define scs_calloc   calloc
#endif

#ifdef DLONG
    #ifdef _WIN64
        typedef __int64 scs_int;
        /* #define scs_int __int64 */
    #else
        typedef long scs_int;
        /* #define scs_int long */
    #endif
#else
    typedef int scs_int;
    /* #define scs_int int */
#endif

#ifndef FLOAT
typedef double scs_float;
#ifndef NAN
#define NAN ((scs_float)0x7ff8000000000000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#else
typedef float scs_float;
#ifndef NAN
#define NAN ((float)0x7fc00000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#endif

#ifdef NULL
#undef NULL
#endif
#define NULL 0

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#ifndef POWF
#ifdef FLOAT
#define POWF powf
#else
#define POWF pow
#endif
#endif

#ifndef SQRTF
#ifdef FLOAT
#define SQRTF sqrtf
#else
#define SQRTF sqrt
#endif
#endif

typedef struct SCS_PROBLEM_DATA Data;
typedef struct SCS_SETTINGS Settings;
typedef struct SCS_SOL_VARS Sol;
typedef struct SCS_INFO Info;
typedef struct SCS_SCALING Scaling;
typedef struct SCS_WORK Work;
typedef struct SCS_CONE Cone;

#endif

