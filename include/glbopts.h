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
#define SCS_VERSION ("3.3.1")

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
#define ACCELERATION_INTERVAL (5)
#define ADAPTIVE_SCALE (1)
#define ADAPTIVE_DIAG_SCALE (1)
#define WRITE_DATA_FILENAME (0)
#define LOG_CSV_FILENAME (0)
#define TIME_LIMIT_SECS (0.)
/* Tolerance to check negativity condition for infeasibility */
#define INFEAS_NEGATIVITY_TOL (1e-9)
/* Number of consecutive residual checks an infeasibility/unboundedness
 * certificate must pass before SCS declares it. Guards against transient
 * iterates (e.g. from acceleration) that momentarily pass a one-shot test. */
#define CERT_PERSISTENCE_CHECKS (2)
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
/* Default AA type: 1 = type-I (better empirical performance, default),
 * 0 = type-II (more numerically stable but typically slower). */
#define ACCELERATION_TYPE_1 (1)
/* Default Tikhonov regularization for the AA least-squares solve. Tuned
 * for type-I; type-II tolerates much smaller (e.g. 1e-12). Users picking
 * type-II will typically lower this. */
#define AA_REGULARIZATION (1e-8)
#define AA_RELAXATION (1.0)
/* Reject AA steps when the output norm exceeds this multiple of the input
 * norm. 1.0 means the AA step must not increase the iterate norm. */
#define AA_SAFEGUARD_FACTOR (1.)
/* Reject AA steps whose weight vector exceeds this norm (prevents
 * numerically unstable extrapolation). */
#define AA_MAX_WEIGHT_NORM (1e10)
/* Max iterative-refinement passes on the γ solve. 0 disables IR; the loop
 * auto-stops once the correction no longer contracts, so this is an upper
 * bound rather than a fixed iteration count. */
#define AA_IR_MAX_STEPS (5)

/* (Dual) Scale updating parameters */
#define MAX_SCALE_VALUE (1e6)
#define MIN_SCALE_VALUE (1e-6)
#define SCALE_NORM NORM /* what norm to use when computing the scale factor */

/* Dynamic diagonal rescaling (stgs->adaptive_diag_scale). Row multipliers
 * move by at most (profile ratio)^DIAG_SCALE_DAMP per update and live in
 * [DIAG_SCALE_MULT_MIN, DIAG_SCALE_MULT_MAX] around the scalar scale. An
 * update fires when the scalar scale updates, or when some damped
 * *clamped* step alone exceeds sqrt(10) (a railed scalar must not freeze
 * the diagonal; a railed multiplier must not keep triggering updates it
 * cannot take). */
#define DIAG_SCALE_DAMP (0.25)
#define DIAG_SCALE_MULT_MIN (1e-3)
#define DIAG_SCALE_MULT_MAX (1e3)
/* Floor on the row-profile denominators, as a fraction of the block's
 * rms denominator (see row_rel_res). Swept over 1e-4..1e-1: every value
 * improves on no floor, 1e-3 is the best on solve count. */
#define DEN_FLOOR_FRAC (1e-3)

/* --- Conjugate gradient (CG) parameters, only used with indirect solver --- */
#define CG_BEST_TOL (1e-12)
/* Each CG solve targets tol = CG_TOL_FACTOR * current_residual. Smaller
 * values give more accurate CG solves at the cost of more CG iterations.
 * With deflation making inner accuracy cheap, 0.03 measured optimal on a
 * netlib / Maros-Meszaros / SOC suite (bracketed on both sides: 0.02
 * starts flipping problems, 0.01 pays more CG for no outer-iteration
 * gain); together with CG_RATE 2.0 it buys roughly 30% fewer outer
 * iterations for roughly 8% more matvecs, which is wall-clock neutral
 * on matvec-dominated problems and favorable on cone-dominated ones.
 * The retuned values target accuracies below the single-precision noise
 * floor, so single-precision builds keep the previous calibration. */
#ifdef SFLOAT
#define CG_TOL_FACTOR (0.2)
#else
#define CG_TOL_FACTOR (0.03)
#endif

/* norm to use when deciding CG convergence */
#ifndef CG_NORM
#define CG_NORM SCS(norm_inf)
#endif
/* cg tol ~ O(1/k^(CG_RATE)); forcing accuracy faster with the iteration
 * count is consumed by the outer loop (Anderson acceleration
 * especially) as fewer iterations. 2.0 is safely interior: 2.25 starts
 * flipping problems. Single precision keeps the previous rate for the
 * same reason as CG_TOL_FACTOR above. */
#ifdef SFLOAT
#define CG_RATE (1.5)
#else
#define CG_RATE (2.0)
#endif
/* Number of approximate small eigenvectors of the preconditioned
 * reduced operator harvested from each cold solve for g = K^{-1}h and
 * used to deflate the warm solves until the next metric change (eigCG,
 * Stathopoulos & Orginos 2010). Costs 2 * DEFLATE_VECTORS * n floats
 * persistent plus a transient Lanczos window during the deep solve;
 * needs USE_LAPACK and is compiled out without it. 0 disables. The
 * environment variables SCS_DEFLATE and SCS_EIGCG_WIN override the
 * count and the window size at runtime. */
#define DEFLATE_VECTORS (30)

#ifdef __cplusplus
}
#endif
#endif
