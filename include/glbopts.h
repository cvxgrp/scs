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
 * the solve, whatever the linear solver, and compares its inf-norm with tol,
 * the smaller inf-norm of SCS's primal and dual residual vectors at the last
 * convergence check, but at least REFINE_TOL_FLOOR. When the residual is
 * above an "on" factor times tol, refinement turns on; once on, it turns off
 * at a later check whose residual is at or below a smaller "off" factor
 * times tol (with one factor for both it switched off and on again, and the
 * unrefined solves in between undid the progress: at 1e-9 that stalled 20
 * of 400 runs with cuDSS). What "on" means depends on the kind of linear
 * solver, chosen at compile time by INDIRECT, so each solver behaves the
 * same in every build and interface:
 *
 * Direct solvers (a factorization: QDLDL, dense, Accelerate, MKL Pardiso,
 * cuDSS). The error of a solve with a factorization is about proportional to
 * the size of the right-hand side. While on, each iteration therefore solves
 * in correction form from the previous iteration's solution u0:
 * K d = rhs - K u0 and u = u0 + d. The iterates move little from one
 * iteration to the next, so rhs - K u0 is far smaller than rhs (10^-3 to
 * 10^-7 of it in the measurements) and u is that much more accurate than a
 * plain solve: on a problem whose plain solves had errors of 1e-7, u stayed
 * at 1e-16 for 70,000 iterations; the errors do not build up, because each
 * solve's error is a contraction of the previous one. The cost is one more
 * multiply by K per iteration, and no extra solve. When the guess is no better
 * than zero (||rhs - K u0|| >= ||rhs||, which happens after a scale update
 * changes K and rhs) that iteration solves plainly instead. At a check while
 * on, a plain solve of the same right-hand side, made only for the decision,
 * decides whether to stay on; the iterate stays the correction-form solution
 * (unless its own residual, one more multiply per check, is clearly larger
 * than the plain solve's: a safeguard for a factorization too inaccurate for
 * the correction form to contract, see REFINE_DIRECT_SAFEGUARD_FACTOR).
 * That matters: a run that put the plain solve in the iterate at every check
 * mixed one inexact solve into every REFINE_CHECK_INTERVAL accurate ones and
 * took up to 31 times the iterations (QFORPLAN at 1e-6 with Accelerate);
 * what slows SCS down is switching between inexact and accurate solves, not
 * the inexact solves themselves. A linear solver can compute the multiplies
 * by K with its own copy of K (scs_lin_sys_residual in linsys.h; cuDSS keeps
 * the vectors on the GPU), else SCS uses the problem data.
 *
 * Indirect solvers (conjugate gradients, CPU and GPU). CG already solves to a
 * tolerance tied to SCS's residuals, and its error is set by that tolerance,
 * not by the size of the right-hand side, so the correction form would not
 * make it more accurate. While on, each iteration's solve is followed by one
 * refinement step: SCS solves K d = rhs - K u (cold, to CG_BEST_TOL) and adds
 * d to u. That costs one more multiply and one more CG solve per iteration;
 * CG fires only when a CG solve misses its tolerance (it stops early or, in
 * single precision, the residual it tracks drifts from the true one).
 *
 * REFINE_CHECK_INTERVAL: the cost while off is one multiply by K and three
 * inf-norms per interval, about 1% of linear-system time with QDLDL (0.4% of
 * the whole solve) and about 2% with cuDSS (0.5% to 8% on the 7 problems
 * measured where it did not fire). 25 is CONVERGED_INTERVAL
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
/* REFINE_DIRECT_ON_FACTOR, REFINE_DIRECT_OFF_FACTOR: direct solvers turn the
 * correction form on when a plain solve's residual is above 0.01 * tol and
 * off when it is at or below 0.001 * tol. Low factors turn it on before the
 * solve's error limits SCS, and keep it on while the residuals fall, so a run
 * rarely switches between inexact and accurate solves. The price is one more
 * multiply per iteration while on, and a plain solve and a multiply at each
 * check while on. Against the previous rule (factors 1 and 0.1, one
 * refinement step instead of the correction form): QDLDL on random LPs and
 * SOCPs at 1e-9 (500 runs) took 2.6% fewer iterations and 2.3% less time;
 * Accelerate on 30 LPs, QPs and SDPs at 1e-4 and 1e-6 (180 runs) took 0.95
 * times the time in geometric mean of the per-run ratios, and 0.89 times
 * the shifted geometric mean of the times (shift 10 s, weighing the long
 * runs more; QFORPLAN at 1e-6: 700 to 975 iterations from 3 starts, instead
 * of 1,300 to 140,900); cuDSS at 1e-9 took
 * 0.86 to 0.91 times the solve time. At the default tolerances (QDLDL, 112
 * problems from 3 starts each) it was on for 1.3% of the iterations, at about
 * 1% more time. Other factors, all in correction form (measured on the
 * prototype): 1 and 0.1 took as many iterations as the previous rule at 1e-9
 * and 0.94 of its time on the Accelerate set (0.89 here), 0.1 and 0.01 1.4%
 * fewer iterations and 0.91; 0.001 and 0.0001 were on for 12.6% of the
 * iterations at the default tolerances, on for most of the solve on more QPs
 * at 1e-8 and 1e-9 (48% of their iterations, against 6.7%), and solved fewer
 * SDPs that stall. A problem that stalls can keep it on for most of its
 * iterations (on 26 that stall at 1e-4 to 1e-6, for most of the solve in 11
 * of 78 runs). In single precision REFINE_TOL_FLOOR puts the off threshold
 * near or below what a single-precision solve reaches, so once on it tends
 * to stay on. */
#define REFINE_DIRECT_ON_FACTOR (0.01)
#define REFINE_DIRECT_OFF_FACTOR (0.001)
/* REFINE_DIRECT_SAFEGUARD_FACTOR, REFINE_DIRECT_SAFEGUARD_FLOOR: at each check
 * while the correction form is on, SCS also computes the correction-form
 * solution's own residual (one more multiply by K per REFINE_CHECK_INTERVAL
 * iterations while on) and replaces the iterate and the next guess with the
 * check's plain solve when that residual is above
 * REFINE_DIRECT_SAFEGUARD_FACTOR times the plain solve's and above
 * REFINE_DIRECT_SAFEGUARD_FLOOR * ||rhs||_inf (or is NaN). This is for a
 * factorization so inaccurate that the correction form does not contract:
 * otherwise the only fallback is ||rhs - K u0|| >= ||rhs||, which lets the
 * error of the iterates grow to the size of rhs first. With both conditions
 * the error at a check is at most about twice a plain solve's.
 * Why not factor 1 and no floor (fall back whenever the correction form is
 * not more accurate): that swapped an equally accurate plain solve into the
 * iterate where both residuals were at rounding level (1.7e-16 of ||rhs||,
 * truss5 at 1e-6 with Accelerate, 9 times in one run) or within 30% of each
 * other (hinf13, 2 times), and each swap perturbs the iterates as the
 * mixing described above does: 3 of 1,240 runs changed (T, TA, A and named
 * QPs), hinf13 from 215,625 to 582,350 iterations. With factor 2 and the
 * floor it fired in none of them, so their iterates are bit for bit those
 * without the safeguard. The floor is 45 (double) to 85 (single) times
 * machine epsilon: a computed residual below it is mostly the rounding of
 * its own computation, so comparing two of them says nothing. */
#define REFINE_DIRECT_SAFEGUARD_FACTOR (2.0)
#ifdef SFLOAT
#define REFINE_DIRECT_SAFEGUARD_FLOOR (1e-5)
#else
#define REFINE_DIRECT_SAFEGUARD_FLOOR (1e-14)
#endif
/* REFINE_INDIRECT_ON_FACTOR: indirect solvers refine when the solve's error
 * is larger than the residuals SCS is driving to zero. While on, each
 * iteration costs one more multiply and one more solve. A smaller factor
 * refines more: 0.1 took slightly fewer iterations at 1e-9 with QDLDL but
 * fired in 15 of the 336 runs at 1e-4 and 1e-6 instead of 9, and where
 * refinement fires it cost more (3% more linear-system time with CG). A
 * larger one leaves solves that limit accuracy unrefined. */
#define REFINE_INDIRECT_ON_FACTOR (1.0)
/* REFINE_INDIRECT_OFF_FACTOR: once on, CG refinement stays on until the
 * solve's error is at or below this fraction of tol. With an off factor of 1
 * the residuals fell to the error of the refined solves, refinement turned
 * off, and until the next check the unrefined solves undid the progress.
 * The price is that refinement stays on longer once it fires: with CG it
 * refined about twice as many solves. In single precision the off threshold
 * is at least 0.1 * REFINE_TOL_FLOOR, near what a refined solve can reach, so
 * this happens more often: at the default tolerances, on 30 LPs and QPs from
 * 3 starts each, refinement was on for most of the solve in 41 of 90 runs
 * with CG (21 with an off factor of 1), mostly runs that stall, where it cost
 * up to 2.1 times the linear-system time per iteration. It also solved more:
 * 47 instead of 41 (gaining 8, losing 2). */
#define REFINE_INDIRECT_OFF_FACTOR (0.1)
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
