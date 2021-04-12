#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

#ifndef SCS
#define SCS(x) scs_##x
#endif

/* SCS VERSION NUMBER ----------------------------------------------    */
#define SCS_VERSION \
  ("2.1.3") /* string literals automatically null-terminated */

/* SCS returns one of the following integers:                           */
#define SCS_INFEASIBLE_INACCURATE (-7)
#define SCS_UNBOUNDED_INACCURATE (-6)
#define SCS_SIGINT (-5)
#define SCS_FAILED (-4)
#define SCS_INDETERMINATE (-3)
#define SCS_INFEASIBLE (-2) /* primal infeasible, dual unbounded   */
#define SCS_UNBOUNDED (-1)  /* primal unbounded, dual infeasible   */
#define SCS_UNFINISHED (0)  /* never returned, used as placeholder */
#define SCS_SOLVED (1)
#define SCS_SOLVED_INACCURATE (2)

/* DEFAULT SOLVER PARAMETERS AND SETTINGS --------------------------    */
#define MAX_ITERS (5000)
#define EPS (1E-5)
#define ALPHA (1.5)
#define RHO_X (1E-3)
#define SCALE (1.0)
#define CG_RATE (2.0)
#define VERBOSE (1)
#define NORMALIZE (1)
#define WARM_START (0)
#define ACCELERATION_LOOKBACK (10)
#define WRITE_DATA_FILENAME (0)

/* redefine printfs and memory allocators as needed */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define scs_printf mexPrintf
#define _scs_free mxFree
#define _scs_malloc mxMalloc
#define _scs_calloc mxCalloc
#define _scs_realloc mxRealloc
#elif defined PYTHON
#include <Python.h>
#include <stdlib.h>
#define scs_printf(...)                              \
  {                                                  \
    PyGILState_STATE gilstate = PyGILState_Ensure(); \
    PySys_WriteStdout(__VA_ARGS__);                  \
    PyGILState_Release(gilstate);                    \
  }
#define _scs_printf printf
#define _scs_free free
#define _scs_malloc malloc
#define _scs_calloc calloc
#define _scs_realloc realloc
#elif (defined(USING_R))
#include <R_ext/Print.h> /* Rprintf etc */
#include <stdio.h>
#include <stdlib.h>
#define scs_printf Rprintf
#define _scs_free free
#define _scs_malloc malloc
#define _scs_calloc calloc
#define _scs_realloc realloc
#else
#include <stdio.h>
#include <stdlib.h>
#define scs_printf printf
#define _scs_free free
#define _scs_malloc malloc
#define _scs_calloc calloc
#define _scs_realloc realloc
#endif

/* Only required for SuiteSparse compatibility: */
#ifndef _scs_printf
#define _scs_printf scs_printf
#endif

#define scs_free(x) \
  _scs_free(x);     \
  x = SCS_NULL
#define scs_malloc(x) _scs_malloc(x)
#define scs_calloc(x, y) _scs_calloc(x, y)
#define scs_realloc(x, y) _scs_realloc(x, y)

#ifdef DLONG
#ifdef _WIN64
#include <stdint.h>
typedef int64_t scs_int;
/* typedef long scs_int; */
#else
typedef long scs_int;
#endif
#else
typedef int scs_int;
#endif

#ifndef SFLOAT
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

#define SCS_NULL 0

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#ifndef POWF
#ifdef SFLOAT
#define POWF powf
#else
#define POWF pow
#endif
#endif

#ifndef SQRTF
#ifdef SFLOAT
#define SQRTF sqrtf
#else
#define SQRTF sqrt
#endif
#endif

#define EPS_TOL (1E-18)
#define SAFEDIV_POS(X, Y) ((Y) < EPS_TOL ? ((X) / EPS_TOL) : (X) / (Y))

#if EXTRA_VERBOSE > 0
#define PRINT_INTERVAL (1)
#define CONVERGED_INTERVAL (1)
#else
/* print summary output every this num iterations */
#define PRINT_INTERVAL (100)
/* check for convergence every this num iterations */
#define CONVERGED_INTERVAL (20)
#endif

/* tolerance at which we declare problem indeterminate */
#define INDETERMINATE_TOL (1e-9)
/* maintain the iterates at this l2 norm (due to homogeneity) */
#define ITERATE_NORM (10.)

#ifdef __cplusplus
}
#endif
#endif
