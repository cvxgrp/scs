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
#define SCS_VERSION ("3.3.2")

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
/* Back-substitution through the KKT factorization loses about eps / rho_x
 * relative accuracy in x, where eps is the machine precision. Single precision
 * cannot afford 1e-6: even min x s.t. x >= 2 stalls at an 8% primal residual.
 * 1e-2 leaves that error around 1e-5, below the default eps_abs and eps_rel. */
#ifdef SFLOAT
#define RHO_X (1E-2)
#else
#define RHO_X (1E-6)
#endif
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
 * update fires when the scalar scale updates, or when some multiplier, as
 * it would be installed (damped, clamped, and after its cone block has been
 * given its form), would move by more than sqrt(10) (a railed scalar must
 * not freeze the diagonal; a railed multiplier must not keep triggering
 * updates it cannot take; a row whose step its block averages away must
 * not trigger one either). */
#define DIAG_SCALE_DAMP (0.25)
#define DIAG_SCALE_MULT_MIN (1e-3)
#define DIAG_SCALE_MULT_MAX (1e3)
/* Floor on the row-profile denominators, as a fraction of the block's
 * rms denominator (see row_rel_res). Swept over 1e-4..1e-1: every value
 * improves on no floor, 1e-3 is the best on solve count. */
#define DEN_FLOOR_FRAC (1e-3)

/* --- Accuracy of the solve for g = K^{-1} [c; -b] ---
 * SCS solves for g at the start of every solve and after every scale update
 * (once per factorization of K) and uses it in every iteration until the
 * next one, so an
 * error in g is repeated in every iteration instead of being averaged out,
 * and it can limit the accuracy SCS reaches. After that solve SCS computes
 * the residual [c; -b] - K g from the problem data, whatever the linear
 * solver, and while the residual is larger than G_SOLVE_TOL times the largest
 * term of the product (and of [c; -b]) it refines g: it solves K d = residual
 * and adds d to g. It stops after G_REFINE_MAX_STEPS steps, or when a step
 * does not at least halve the residual, and discards a step that makes it
 * larger. Tradeoff: one multiply by K per solve for g, plus one solve and one
 * multiply per refinement step. For the direct solvers that is small next to
 * the factorization; when the solve is already that accurate only the
 * multiply is paid. CG stops on an absolute residual of its reduced system
 * (CG_BEST_TOL), which on badly scaled data can leave a much larger relative
 * residual here, so CG's g is checked and refined too, and for CG a
 * refinement step is another cold CG solve, which can cost up to about as
 * many iterations as the solve for g itself. Single precision cannot reach
 * 1e-12, so it uses 1e-5. */
#ifdef SFLOAT
#define G_SOLVE_TOL (1e-5)
#else
#define G_SOLVE_TOL (1e-12)
#endif
#define G_REFINE_MAX_STEPS (3)

/* --- Accuracy of each iteration's solve ---
 * Every iteration solves K u = [R_x v_x; -R_y v_y] (K as above), and SCS's
 * residuals cannot fall much below the error of that solve. On iterations
 * that are multiples of REFINE_CHECK_INTERVAL SCS computes the residual of
 * the solve from the problem data, whatever the linear solver, and compares
 * its inf-norm with tol, the smaller inf-norm of SCS's primal and dual
 * residual vectors at the last convergence check, but at least
 * REFINE_TOL_FLOOR. If it is larger than REFINE_FACTOR * tol, SCS solves
 * K d = residual and adds d to u, and then checks every iteration until a
 * residual is at or below REFINE_OFF_FACTOR * tol. In double precision at
 * default tolerances this almost never fires (in 9 of 336 benchmark runs at
 * 1e-4 and 1e-6, on 0.01% of solves); at 1e-8 or 1e-9 a direct solver's
 * error can be as large as the residuals, and refining is what lets SCS
 * converge. CG
 * already solves to a fraction of tol, so with CG it fires only when a CG
 * solve misses that (it stops early or, in single precision, the residual it
 * tracks drifts from the true one).
 *
 * REFINE_CHECK_INTERVAL: the cost while off is one multiply by K and three
 * inf-norms per interval, about 1% of linear-system time with QDLDL (0.4% of
 * the whole solve) and about 2% with cuDSS. 25 is CONVERGED_INTERVAL
 * (except in debug builds with VERBOSITY > 0), so each check has a tol the
 * previous one did not. That tol
 * is 25 iterations old: the solve of iteration i comes before the
 * convergence check of iteration i, so the check at iteration 25k compares
 * with the residuals of iteration 25(k-1). A shorter interval catches an
 * inaccurate solve sooner and costs proportionally more. 50 halved the cost
 * but caught inaccurate solves later: at 1e-9 it took more iterations than 25
 * in 49 of 500 runs with QDLDL (fewer in 4) and in 134 of 400 with cuDSS
 * (fewer in 53). */
#define REFINE_CHECK_INTERVAL (25)
/* REFINE_FACTOR: refine when the solve's error is larger than the residuals
 * SCS is driving to zero. While on, each iteration costs one more multiply
 * and one more solve. A smaller factor refines more: 0.1 took slightly fewer
 * iterations at 1e-9 but fired in 15 of the 336 runs at 1e-4 and 1e-6
 * instead of 9, and where refinement fires it cost more (4% more
 * linear-system time on stalled SDPs, 3% with CG). A larger one leaves
 * solves that limit accuracy unrefined. The alternatives measured worse:
 * refining every solve costs 2.2x linear-system time and lost problems, and
 * refining near convergence refined 17% of solves at default tolerances.
 * Single precision uses the same factor: there a QDLDL solve's error is
 * about the size of the residuals at the default tolerances, and refining is
 * what lets SCS converge (test_solve_refine_tight). */
#define REFINE_FACTOR (1.0)
/* REFINE_OFF_FACTOR: once on, refinement stays on until the solve's error is
 * at or below this fraction of tol. Turning it off as soon as the error is
 * back below tol (an off factor of REFINE_FACTOR) made it flicker: the
 * residuals fall to the error of the refined solves, refinement turns off,
 * and until the next check up to REFINE_CHECK_INTERVAL iterations later the
 * unrefined solves undo the progress. At 1e-9 that stalled 20 of 400 runs at
 * max_iters with cuDSS, and with QDLDL 12 of 500 runs took at least twice the
 * iterations they took with 0.1 (up to 34 times). With 0.1 no cuDSS run
 * stalled (1,000 runs, with and without Anderson acceleration) and none of
 * the 500 QDLDL runs took more iterations. On other problems it moved
 * chaotic runs both ways: with QDLDL at 1e-8 and 1e-9, 0.1 took more
 * iterations in 4 of 252 runs (up to 36% more) and fewer in 14, and on SDPs
 * that stall at 1e-6 more in 9 of 62 (up to 3.3 times) and fewer in 3. The
 * price is that refinement stays on longer once it fires: with CG it refined
 * about twice as many solves, and on a problem that stalls it can stay on
 * for most of the remaining iterations (an SDP that ran to max_iters at 1e-6
 * refined 64% of its solves), at up to one more solve and multiply per
 * iteration. In single precision the off threshold is at least 0.1 *
 * REFINE_TOL_FLOOR, near what a refined solve can reach, so this happens
 * more often: at the default tolerances, on 30 LPs and QPs from 3 starts
 * each, refinement was on for most of the solve in 27 of 90 runs with QDLDL
 * (12 with an off factor of 1) and in 41 with CG (21), mostly runs that
 * stall, where it cost up to 2.1 times the linear-system time per iteration.
 * It also solved more: 60 instead of 57 with QDLDL, and 47 instead of 41
 * with CG (gaining 8, losing 2). A smaller value keeps it on
 * longer still; values between 0.1 and 1 were not measured. */
#define REFINE_OFF_FACTOR (0.1)
/* REFINE_TOL_FLOOR: the smallest tol. Refinement cannot take a residual
 * much below rounding level, so below this floor it would chase rounding
 * errors on every remaining iteration: when one of SCS's residual vectors is
 * zero (an infeasibility or unboundedness certificate: without the floor it
 * fired on test-suite certificates with every linear solver, and with CG
 * stayed on for 75 of 100 iterations of test_unbounded_inaccurate), or in
 * single precision when the residuals near rounding level. A higher floor
 * skips refinement that could still help, a lower one refines in vain.
 * Double precision uses 1e-12 (CG_BEST_TOL, the accuracy of SCS's most
 * accurate solves), which changed none of the 520 benchmark runs it was
 * checked on. Single precision uses 1e-5: the single-precision solves that
 * need refinement have tol near 1e-3, and 1e-6 refined more without
 * converging faster (test_adaptive_scale with CG, measured with an off
 * factor of 1: 54 refinement steps and 225 iterations, against 22 and 175
 * with 1e-5; that test takes 175 to 275 iterations from nearby initial
 * scales). */
#ifdef SFLOAT
#define REFINE_TOL_FLOOR (1e-5)
#else
#define REFINE_TOL_FLOOR (1e-12)
#endif

/* --- Conjugate gradient (CG) parameters, only used with indirect solver --- */
/* Accuracy of the most accurate CG solves (the cold solve for g, and the floor
 * of the warm ones). Single precision cannot reach 1e-12: a CG asked for it
 * runs to its iteration limit on every solve and its rounding errors then
 * reach the iterates (issue_499's infeasibility certificate stalled with
 * tau at rounding level). 1e-6 is about ten times float's machine epsilon;
 * 1e-7 also passes the test suite, 1e-5 starts losing accuracy. */
#ifdef SFLOAT
#define CG_BEST_TOL (1e-6)
#else
#define CG_BEST_TOL (1e-12)
#endif
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
