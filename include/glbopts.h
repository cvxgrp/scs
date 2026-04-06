/*
 * Global options, default parameter values, and platform-specific macros.
 *
 * Defines printing/memory allocation macros (adapts to MATLAB, Python, R),
 * math precision macros (float vs double), default solver constants, and
 * internal algorithm tuning parameters. This is an internal header; the
 * public API is in scs.h.
 */

#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include <math.h>

#ifndef SCS
#define SCS(x) _scs_##x
#endif

/* SCS VERSION NUMBER ----------------------------------------------     */
/* string literals automatically null-terminated */
#define SCS_VERSION ("3.2.11")

/* verbosity level */
#ifndef VERBOSITY
#define VERBOSITY (0)
#endif

/* DEFAULT SOLVER PARAMETERS AND SETTINGS --------------------------     */
/* If you update any of these you must update the documentation manually */
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
/* Tolerance to check negativity condition for infeasibility */
#define INFEAS_NEGATIVITY_TOL (1e-9)
/* redefine printfs as needed */
#if NO_PRINTING > 0     /* Disable all printing */
#define scs_printf(...) /* No-op */
#else
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define scs_printf mexPrintf
#elif defined PYTHON
#include <Python.h>
/* see:
 * https://cython-users.narkive.com/jRjjs3sK/reacquire-gil-for-printing-in-wrapped-c-library
 */
#define scs_printf(...)                                                        \
  {                                                                            \
    PyGILState_STATE gilstate = PyGILState_Ensure();                           \
    PySys_WriteStdout(__VA_ARGS__);                                            \
    PyGILState_Release(gilstate);                                              \
  }
#elif defined R_LANG
#include <R_ext/Print.h> /* Rprintf etc */
#include <stdio.h>
#include <stdlib.h>
#define scs_printf Rprintf
#else
#include <stdio.h>
#include <stdlib.h>
#define scs_printf printf
#endif
#endif

/* redefine memory allocators as needed */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
/* Use mexMakeMemoryPersistent so allocations survive across MEX calls.
 * Required for the workspace API (scs_init/scs_solve/scs_finish). */
static inline void *_scs_mex_malloc(size_t n) {
  void *p = mxMalloc(n);
  if (p) mexMakeMemoryPersistent(p);
  return p;
}
static inline void *_scs_mex_calloc(size_t count, size_t size) {
  void *p = mxCalloc(count, size);
  if (p) mexMakeMemoryPersistent(p);
  return p;
}
static inline void *_scs_mex_realloc(void *ptr, size_t n) {
  void *p = mxRealloc(ptr, n);
  if (p) mexMakeMemoryPersistent(p);
  return p;
}
#define scs_free mxFree
#define scs_malloc _scs_mex_malloc
#define scs_calloc _scs_mex_calloc
#define scs_realloc _scs_mex_realloc
#elif defined PYTHON
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define scs_free PyMem_RawFree
#define scs_malloc PyMem_RawMalloc
#define scs_realloc PyMem_RawRealloc
#define scs_calloc PyMem_RawCalloc
#else
#define scs_free PyMem_Free
#define scs_malloc PyMem_Malloc
#define scs_realloc PyMem_Realloc
static inline void *scs_calloc(size_t count, size_t size) {
  void *obj = PyMem_Malloc(count * size);
  memset(obj, 0, count * size);
  return obj;
}
#endif
#elif defined R_LANG
#include <stdio.h>
#include <stdlib.h>
#define scs_free free
#define scs_malloc malloc
#define scs_calloc calloc
#define scs_realloc realloc
#else
#include <stdio.h>
#include <stdlib.h>
#define scs_free free
#define scs_malloc malloc
#define scs_calloc calloc
#define scs_realloc realloc
#endif

#ifndef SFLOAT
#ifndef NAN
#define NAN (HUGE_VAL - HUGE_VAL)
#endif
#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif
#else
#ifndef NAN
#define NAN ((float)(HUGE_VAL - HUGE_VAL))
#endif
#ifndef INFINITY
#define INFINITY ((float)HUGE_VAL)
#endif
#endif

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifdef SFLOAT
#define SQRTF sqrtf
#define EXPF expf
#define LOGF logf
#define ABS fabsf
#define POWF powf
#else
#define SQRTF sqrt
#define EXPF exp
#define LOGF log
#define ABS fabs
#define POWF pow
#endif

#ifdef DLONG
#define IABS llabs
#else
#define IABS abs
#endif

/* Force SCS to treat the problem as (non-homogeneous) feasible for this many
 * iterations. Acts like an implicit warm-start biased towards feasibility,
 * which is the most common use-case. During these iterations tau is fixed
 * at 1 and kappa is fixed at 0. */
#define FEASIBLE_ITERS (1)

/* Minimum iterations between heuristic scale updates. Prevents scale
 * from changing too frequently before the iterates have stabilized. */
#define RESCALING_MIN_ITERS (100)

#define _DIV_EPS_TOL (1E-18)
#define SAFEDIV_POS(X, Y)                                                      \
  ((Y) < _DIV_EPS_TOL ? ((X) / _DIV_EPS_TOL) : (X) / (Y))

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

/* Factor which scales the tau diagonal entry in the linear system.
 * Larger values stabilize tau but slow convergence. 10 is a good balance
 * for most problems. */
#define TAU_FACTOR (10.)

/* --- Anderson Acceleration (AA) parameters --- */
#define AA_RELAXATION (1.0)
#define AA_REGULARIZATION_TYPE_1 (1e-6)
#define AA_REGULARIZATION_TYPE_2 (1e-10)
/* Reject AA steps when the output norm exceeds this multiple of the input
 * norm. 1.0 means the AA step must not increase the iterate norm. */
#define AA_SAFEGUARD_FACTOR (1.)
/* Reject AA steps whose weight vector exceeds this norm (prevents
 * numerically unstable extrapolation). */
#define AA_MAX_WEIGHT_NORM (1e10)

/* (Dual) Scale updating parameters */
#define MAX_SCALE_VALUE (1e6)
#define MIN_SCALE_VALUE (1e-6)
#define SCALE_NORM NORM /* what norm to use when computing the scale factor */

/* --- Conjugate gradient (CG) parameters, only used with indirect solver --- */
#define CG_BEST_TOL (1e-12)
/* Each CG solve targets tol = CG_TOL_FACTOR * current_residual. Smaller
 * values give more accurate CG solves at the cost of more CG iterations. */
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
