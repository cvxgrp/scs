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
#define SCS_VERSION                                                            \
  ("3.0.0") /* string literals automatically null-terminated */

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

/* verbosity level */
#ifndef VERBOSITY
#define VERBOSITY (0)
#endif

/* DEFAULT SOLVER PARAMETERS AND SETTINGS --------------------------    */
#define MAX_ITERS (100000)
#define EPS_REL (1E-4)
#define EPS_ABS (1E-4)
#define EPS_INFEAS (1E-7)
#define ALPHA (1.5)
#define RHO_X (1E-6)
#define SCALE (0.1)
#define VERBOSE (1)
#define NORMALIZE (1)
#define WARM_START (0)
#define ACCELERATION_LOOKBACK (10)
#define ACCELERATION_INTERVAL (10)
#define ADAPTIVE_SCALE (1)
#define WRITE_DATA_FILENAME (0)
#define LOG_CSV_FILENAME (0)
#define TIME_LIMIT_SECS (0.)

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
#define scs_printf(...)                                                        \
  {                                                                            \
    PyGILState_STATE gilstate = PyGILState_Ensure();                           \
    PySys_WriteStdout(__VA_ARGS__);                                            \
    PyGILState_Release(gilstate);                                              \
  }
/* only for SuiteSparse */
#define _scs_printf PySys_WriteStdout
#if PY_MAJOR_VERSION >= 3
#define _scs_free PyMem_RawFree
#define _scs_malloc PyMem_RawMalloc
#define _scs_realloc PyMem_RawRealloc
#define _scs_calloc PyMem_RawCalloc
#else
#define _scs_free PyMem_Free
#define _scs_malloc PyMem_Malloc
#define _scs_realloc PyMem_Realloc
static inline void *_scs_calloc(size_t count, size_t size) {
  void *obj = PyMem_Malloc(count * size);
  memset(obj, 0, count * size);
  return obj;
}
#endif
#elif defined R_LANG
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

#define scs_free(x)                                                            \
  _scs_free(x);                                                                \
  x = SCS_NULL
#define scs_malloc(x) _scs_malloc(x)
#define scs_calloc(x, y) _scs_calloc(x, y)
#define scs_realloc(x, y) _scs_realloc(x, y)

#ifdef DLONG
/*#ifdef _WIN64
#include <stdint.h>
typedef int64_t scs_int;
#else
typedef long scs_int;
#endif
*/
typedef long long scs_int;
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

/* Force SCS to treat the problem as (non-homogeneous) feasible for this many */
/* iters. This acts like a warm-start that biases towards feasibility, which */
/* is the most common use-case */
#define FEASIBLE_ITERS (1)

/* how many iterations between heuristic residual rescaling */
#define RESCALING_MIN_ITERS (100)

#define EPS_TOL (1E-18)
#define SAFEDIV_POS(X, Y) ((Y) < EPS_TOL ? ((X) / EPS_TOL) : (X) / (Y))

#if VERBOSITY > 0
#define PRINT_INTERVAL (1)
#define CONVERGED_INTERVAL (1)
#else

/* print summary output every this num iterations */
#define PRINT_INTERVAL (250)
/* check for convergence every this num iterations */
#define CONVERGED_INTERVAL (25)
#endif

/* maintain the iterates at L2 norm =  ITERATE_NORM * sqrt(n+m+1) */
#define ITERATE_NORM (1.)

/* Which norm to use for termination checking etc */
/* #define NORM SCS(norm_2) */
#define NORM SCS(norm_inf)

/* Factor which is scales tau in the linear system update */
/* Larger factors prevent tau from moving as much */
#define TAU_FACTOR (10.)

/* Anderson acceleration parameters: */
#define AA_RELAXATION (1.0)
#define AA_REGULARIZATION_TYPE_1 (1e-6)
#define AA_REGULARIZATION_TYPE_2 (1e-10)
/* Safeguarding norm factor at which we reject AA steps */
#define AA_SAFEGUARD_FACTOR (1.)
/* Max allowable AA weight norm */
#define AA_MAX_WEIGHT_NORM (1e10)

/* (Dual) Scale updating parameters */
#define MAX_SCALE_VALUE (1e6)
#define MIN_SCALE_VALUE (1e-6)
#define SCALE_NORM NORM /* what norm to use when computing the scale factor */

/* CG == Conjugate gradient */
/* Linear system tolerances, only used with indirect */
#define CG_BEST_TOL (1e-12)
/* This scales the current residuals to get the tolerance we solve the
 * linear system to at each iteration. Lower factors require more CG steps
 * but give better accuracy */
#define CG_TOL_FACTOR (0.2)

/* norm to use when deciding CG convergence */
#ifndef CG_NORM
#define CG_NORM SCS(norm_inf)
#endif
/* cg tol ~ O(1/k^(CG_RATE)) */
#define CG_RATE (1.5)

#ifdef __cplusplus
}
#endif
#endif
