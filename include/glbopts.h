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
#include <float.h>
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

/* The R_y entry (rho_y) of an equality row, a zero-cone row, is
 * 1 / (ZERO_CONE_RY_FACTOR * scale); a row of any other cone gets 1 / scale
 * (see SCS(set_r_y) in cones.c). The zero cone's dual is free, so a small
 * rho_y there costs nothing in the cone projection and lets the linear system
 * solve determine those rows almost entirely. Tradeoff: the larger the factor,
 * the more tightly the equalities are held each iteration, and the more the
 * linear system, rather than the projection, has to move them. The adaptive
 * per-row multipliers (DIAG_SCALE_MULT_MIN/MAX) act around this base. */
#define ZERO_CONE_RY_FACTOR (1000.)

/* BOX_T_RY_FACTOR: the weight of the box cone's first entry t. It is an
 * equality row's, ZERO_CONE_RY_FACTOR (1000), in single and double precision.
 *
 * What it is. The R_y entry (rho_y) of t is 1 / (BOX_T_RY_FACTOR * scale),
 * where a row of any cone other than the zero cone gets 1 / scale; with
 * adaptive_diag_scale the row's own multiplier divides it further, as it
 * does every row's (see SCS(set_r_y) in cones.c). It applies to box cones of
 * size 2 or more. A box cone of size 1 is t >= 0 alone, a nonnegative entry
 * with no bounds, and keeps 1 / scale.
 *
 * Why an equality row's weight. In the box cone's usual form, l <= a'x <= u
 * is written as (1, b - A x) in {(t, s) : t l <= s <= t u}: the row of A for
 * t is zero and its b entry is 1, so the row says s_t = tau, an equality,
 * although t is a cone entry. Weighted like an inequality row (1 / scale, as
 * in SCS 3.3), t could stop several percent away from 1 while SCS's residual
 * still passed: that residual weighs t - 1 against the largest entry of A x,
 * but every bound is met only up to |t - 1| times the bound, so problems
 * with large bounds stopped with their bounds violated. Weighted like an
 * equality row, t is held to tau as an equality is.
 *
 * Where it acts (SCS(set_r_y) has the algebra):
 * - The linear system. A_t = 0, so the t row is decoupled from x and gives
 *   y_t = v_t - tau b_t / rho_y.
 * - tau's update. The t entry adds b_t^2 / rho_y = b_t^2 * factor * scale
 *   (times the row's multiplier) to the leading coefficient of the quadratic
 *   that root_plus solves for tau, against tau's own weight TAU_FACTOR = 10.
 *   Where that term dominates, tau b_t is anchored to rho_y v_t, the scaled
 *   t entry of v, which moves by O(rho_y) per iteration: the weight pins tau.
 * - The box projection. Its metric weighs t by 1 / rho_y, so the larger the
 *   factor, the nearer to tau b_t the t at which the box rows are clipped.
 *
 * Tradeoff. Larger: t tracks tau more tightly, and fewer stops claim solved
 * or infeasible at points that fail an independent check (the table below).
 * But the pin on tau grows with the factor, so tau moves more slowly:
 * - Certificates of infeasibility take more iterations. Against 1000, the
 *   infeasible problems below took 1.04x the iterations at 100, 1.12x at 1e4
 *   and 1.69x at 1e5. As the factor goes to infinity tau freezes and SCS
 *   becomes Douglas-Rachford with tau fixed: it still converges when a
 *   solution exists, but certificates come only at rate O(1/k).
 * - A problem that needs tau to move can stall. QPLIB 8515 as a QP benchmark
 *   harness poses it, at eps 1e-4 with QDLDL: its normalized b_t is 9.46, so
 *   at the default scale the t row's pull on tau is 8.9e3 at 1000, against
 *   TAU_FACTOR = 10 (8.9 at factor 1). With tau pinned the dual objective
 *   lags, the first convergence check (iteration 25) misses on the gap, and
 *   then Anderson acceleration extrapolates a slow primal drift and the
 *   adaptive scale walks to its floor. From the default start every factor
 *   from 10 to 1e5 reaches a 20,000-iteration cap. At 1000 so does every
 *   start within 4 ulps of the default scale (#548's code solved 2 of those
 *   9, in 13,925 and 16,950 iterations), and the harness's 300 s limit stops
 *   it at 685,000 iterations with the gap at 5.6e-2. Only factor 1, the
 *   weight in SCS 3.3, solves it, in 25 iterations as SCS 3.3.2 does. 3.3.2
 *   itself reaches the 20,000 cap at eps 2e-5, 1e-5 and 1e-6 (it solves at
 *   3e-5): it solves 8515 only by that early stop.
 * Smaller: t drifts from tau, and fewer problems solve (factor 1 below). Any
 * positive value is a valid metric: it changes the path SCS takes, never the
 * solution.
 *
 * Data. 492 benchmark problems, 483 of them with box cones (Netlib,
 * Kennington, Maros-Meszaros, Netlib infeasible and MIPLIB LP relaxations),
 * each at eps 1e-4 and 1e-6: 984 runs per factor, QDLDL, one thread, every
 * status checked independently against the problem's KKT conditions at
 * 10 x eps. The problems were split into two halves, one to choose and one
 * to confirm, before any run. One run each at the default start:
 *
 *  factor  solved  pen. SGM [95%]        >= 1 s  rejected  infeas. it  |t-1|
 *       1     922  1.093 [0.968, 1.225]   0.995        30         465  1.2e-3
 *      10     930  1.000 [0.895, 1.111]   0.925        23         381  2.1e-4
 *     100     938  0.936 [0.848, 1.032]   0.897        19         429  7.4e-5
 *    1000     933  1                      1            17         412  1.2e-5
 *     1e4     940  0.910 [0.827, 0.991]   0.960        12         462  5.0e-6
 *     1e5     943  0.874 [0.769, 0.974]   0.996        12         698  1.3e-6
 *
 * Then factors 100 to 1e5 ran 8 more starts (the initial scale moved by +-1
 * to +-4 ulps) on the 163 problem / eps pairs where, at the default start,
 * 100, 1e4 or 1e5 had differed from 1000 in verified outcome or in total
 * wall time by more than 2x either way; each pair counts once, as its mean
 * over its 9 starts. These shrink the differences above, as expected for
 * pairs picked because they differed:
 *
 *   factor  expected  pen. SGM [95%]        >= 1 s  rejected
 *             solved
 *      100    936.78  0.984 [0.926, 1.048]   0.946       142
 *     1000    936.22  1                      1           124
 *      1e4    938.78  0.948 [0.890, 1.003]   1.003        87
 *      1e5    941.78  0.914 [0.828, 0.993]   1.029        78
 *
 * - solved: verified solved or certified infeasible, of 984; expected solved
 *   counts each of the 163 pairs as the fraction of its starts solved.
 * - pen. SGM: shifted geometric mean of the time (shift 10 s, a failure
 *   charged 3 x the time limit), over factor 1000's; the interval is a
 *   bootstrap over problems.
 * - >= 1 s: geometric mean of the time over factor 1000's, on the runs both
 *   solve whose geometric-mean time is 1 s or more.
 * - rejected: runs that stopped claiming solved or infeasible where the
 *   independent check rejected the claim (second table: over all starts).
 * - infeas. it: shifted geometric mean (shift 10) of the iterations on the
 *   58 runs of the 29 Netlib infeasible problems; 51 certificates verify at
 *   every factor from 100 up, 50 at 1 and 10.
 * - |t-1|: its 90th percentile where SCS stops, at eps 1e-4, by nearest rank
 *   (the ceil(0.9 n)-th smallest), over the n = 453 to 455 runs with a box
 *   cone of size 2 or more that return a t (a stop that claims
 *   infeasibility returns none).
 *
 * Noise. In verified solves, 100, 1e4 and 1e5 are all within noise of 1000.
 * At the default start the largest difference, 1e5's +10, has sign test
 * p = 0.076; over the ensembles every bootstrap interval for the difference
 * includes 0, over all problems and on each half (100: +0.56 [-4.78, +6.00],
 * 1e4: +2.56 [-2.22, +7.89], 1e5: +5.56 [-2.22, +14.11]). Factor 10 is
 * within noise too (3 fewer). Factor 1 is worse beyond it: 10 fewer solved
 * than 1000 on the confirming half (sign test p = 0.031; 11 fewer over all,
 * p = 0.061), and 16 to 21 fewer than 100, 1e4 and 1e5 (p <= 0.0025).
 *
 * Why 1000. With 100 to 1e5 within noise, the choice is holistic. 1000 is
 * the weight of an equality row, which the t row is, rather than a value
 * picked as the best of six on the same data: the ensembles shrank the
 * default-start gains of 100, 1e4 and 1e5 over 1000 (+5, +7 and +10 solved)
 * to +0.56, +2.56 and +5.56, so the leaders' edges on this data are likely
 * overstated. It is the same value in single and double precision. It is
 * the value #548 was tested and benchmarked with. And it pins tau 10 times
 * less than 1e4. Against 1e4 it gives up 2.56 expected solves (within
 * noise), with a penalized SGM 1 / 0.948 = 1.05 times 1e4's and 124 rejected
 * claims against 87, while long runs are level (1e4: 1.003) and 1e4 takes
 * 1.12 times its iterations on the infeasible problems. 1000 was the
 * maintainer's holistic choice within noise.
 *
 * Known cases at 1000, solved of the 9 starts (1e4 and 1e5 in brackets):
 * - QGROW22: 0 at each eps (2, 9). Every start stops at 250 iterations
 *   claiming dual infeasibility, with t collapsed to about 1e-8, and the
 *   check rejects the claim.
 * - QGROW15 at eps 1e-4: 0 (2, 9); on 8 starts t ends between 1.7e3 and
 *   1.4e5. QGROW7 at 1e-4: 4 (9, 9).
 * - KSIP at 1e-6: 0, all at the iteration cap (1, 0; factor 100: 7).
 * - QGFRDXPN: 3 at each eps, the others at the iteration cap (9, 9; 100: 9).
 * - pilotnov at 1e-6: 0, all at the iteration cap (0, 9).
 * - QPLIB 8515 at 1e-4 (above).
 * Where 1000 does better: mzzv42z at 1e-4, 8 (1, 1), and mzzv11 at 1e-4, 9
 * (3, 4), nearly all accepted by the check at the time limit; QSEBA at
 * 1e-6, 9 (7, 6).
 *
 * Single precision uses the same value. Before, it kept 1000 while double
 * precision had 1e4, because 1e4 failed there: it makes tau 10 times
 * smaller at the solution, and the absolute floor on CG's tolerance
 * (CG_BEST_TOL, 1e-6 in single precision) then kept CG from improving the
 * iterate. hs21_tiny_qp's last solve as the test runs it (CG, eps 1e-4,
 * normalize and adaptive_scale off, USE_LAPACK=0) ended with tau 0.0026
 * against 0.026 at 1000, and reached max_iters from 23 of 41 starting
 * scales near the default on x86 (-march=x86-64-v3, AMD and Intel) and 26
 * of 41 on an Apple M4, so CI's single-precision tests failed on AMD
 * runners. With QDLDL at eps 1e-5 and default settings it stalled from 8
 * (AMD) and 17 (M4) of 41. At 1000, none of either. A CG floor of 1e-7
 * removed the CG stall at 1e4, while at 1000 even 1e-5 causes none. The
 * absolute floor remains a weakness: any change that makes tau smaller in
 * single precision can stall CG again.
 *
 * What it needs. The box projection must be exact (proj_box_cone): with the
 * weight, a projection that stopped short of t* made ADMM cycle. The weight
 * also makes tau smaller at the solution, since the t entry of v is
 * y_t + s_t / rho_y with s_t = tau b_t (hs21_tiny_qp's first solve, in double
 * precision with the test's settings and acceleration: tau 0.019 at 1000,
 * 0.34 at factor 1). That is harmless in double precision, but fixed
 * absolute floors on the iterate, such as the CG floor in single precision,
 * count for more (docs/src/algorithm/scale.rst).
 *
 * Not measured: the GPU and cuDSS solvers, MKL, timings on x86 (the sweep
 * and the ensembles ran on 64-bit ARM Linux), single precision beyond the
 * cases above, and problems outside the benchmark set above other than QPLIB
 * 8515. The sweep and the ensembles ran on SCS #548 with only this value
 * changed. Since then #549 refined the solve for g further (G_TOL_EPS),
 * which moves g in double precision only at rounding level, and nothing
 * that touches the box cone has changed. */
#define BOX_T_RY_FACTOR (ZERO_CONE_RY_FACTOR)

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

/* SCS_MACH_EPS: the machine epsilon of scs_float (DBL_EPSILON, or FLT_EPSILON
 * with SFLOAT), the unit of the tolerance of the solve for g below and of
 * REFINE_ROUNDING_EPS. */
#ifdef SFLOAT
#define SCS_MACH_EPS (FLT_EPSILON)
#else
#define SCS_MACH_EPS (DBL_EPSILON)
#endif

/* --- Accuracy of the solve for g = K^{-1} [c; -b] ---
 * SCS solves for g at the start of every solve and after every scale update
 * (once per factorization of K) and uses it in every iteration until the
 * next one, so an error in g is repeated in every iteration instead of being
 * averaged out, and it can limit the accuracy SCS reaches. g is therefore
 * solved to about machine precision, with the linear solve's own warm start
 * and refinement (scs_solve_lin_sys, linsys.h; update_work_cache in
 * src/scs.c). A plain solve (CG: from zero to CG_BEST_TOL) gives g0 and its
 * residual; if that is above
 *   tol_g = G_TOL_EPS * SCS_MACH_EPS * || |[c; -b]| + |K| |g0| ||
 * (absolute values entry by entry: each row's sum of the absolute values of
 * the terms that row of the residual is computed from, and the largest of
 * these over the rows), the solve is repeated from g0 to tol_g, with at
 * most G_REFINE_MAX_STEPS steps in all: for a direct solver
 * the correction form from g0 is the first, for CG a CG solve from g0 to
 * tol_g; each further step solves for the residual, computed with K, and
 * adds the correction (CG: by CG from zero to tol_g). A step is kept only if
 * it lowers the residual, the steps stop after one that does not halve it
 * (LIN_SYS_REFINE_CONTRACTION), and g0 is kept unless the result is better.
 * Cost: one multiply by K per solve for g (the residual of g0) and a pass
 * over A and P for the size of the terms; when g0 misses tol_g, a solve and
 * a multiply per step, and for a direct solver one multiply more (the
 * correction form computes the residual of g0 again). For CG a step costs
 * the CG iterations that reach tol_g: in all, about 14% of the first solves'
 * on random SOCPs in double precision (60% in single), but where CG cannot
 * reach tol_g each step runs to CG's iteration limit (up to five times the
 * first solve, on badly scaled QPs).
 *
 * G_TOL_EPS: just below the residual that refinement reaches, so that the steps
 * end at the halving stop, where the linear solver can do no better, rather
 * than at tol_g: g to about machine precision, modulo what the solver achieves.
 * Refinement run to that stop ends at 0.3-0.4 machine epsilons of the size of
 * the terms (median; 99th percentile 0.6-1.1, at most 1.7, over about 30,000
 * solves in the A/B of this constant and its diagnosis), with QDLDL and CPU
 * and GPU CG in double and single precision (random SOCPs, residuals computed
 * in double): the rounding errors of a row's terms partly cancel, so the
 * computed residual rounds by far less than an epsilon per term. (The "few
 * epsilons" of rounding level elsewhere, REFINE_ROUNDING_EPS, kkt_terms_inf
 * in src/scs.c and linsys.h, is a bound with room; this is the typical size.)
 * A tol_g the solver can meet stops CG early, because a CG step solves only to
 * tol_g and lands just below it. In single precision that costs time and
 * solves: an iteration that uses a computed g is exact for the data K g, that
 * is for c - r_x and b + r_y (r = [c; -b] - K g), so g's residual is a floor
 * under SCS's residuals, and at eps 1e-6 (about eight single-precision
 * epsilons) it can be above SCS's tolerance where rows are dense. With
 * G_TOL_EPS 4, single-precision CG's g stopped at 3.4 epsilons of its terms
 * (median) on random SOCPs with 100 variables and dense rows, and 72 of 100 ran
 * to max_iters or the time limit; with 0.25, none did, and over all 500 random
 * SOCPs single-precision CG was about three times faster
 * (test_g_accuracy_iteration solves one like them). With 1, CG's g stopped at
 * 0.8 epsilons and some runs still reached the floor; lower than 0.25 gives the
 * same g with more rejected steps. In double precision the floor is about 1e-6
 * of SCS's tolerance at eps 1e-9 either way. Cost: usually one step more per
 * solve for g, the one that does not halve the residual (a solve and a
 * multiply); with CG in double precision, up to 2% more CG iterations and
 * 1-1.5% more time in all on random SOCPs (QDLDL: 0.5%). Where CG's steps run
 * to its iteration limit the extra step costs most: on badly scaled QPs in
 * single precision (HUES-MOD and QSCFXM1, 9 starts each) the steps took 2.2
 * times the CG iterations they took with G_TOL_EPS 4, and the runs 26% more
 * CG iterations in all, while 15 of the 18 runs solved (KKT-verified)
 * instead of 9.
 *
 * The size is taken row by row, not as ||[c; -b]|| + ||K|| ||g0||
 * (the normwise backward error, which pairs K's largest row with g0's
 * largest entry even where they do not meet): on badly scaled data that was
 * up to 1e9 times larger (regress340, QDLDL), and let most solves for g stop
 * with residuals up to 1e5 machine epsilons of their terms. It is one
 * number, the largest row's, because the solve's tol bounds the largest
 * entry of the residual: a row whose terms are much smaller than the
 * largest row's is held to the largest row's rounding, not to its own.
 *
 * G_REFINE_MAX_STEPS: LAPACK's ITMAX (xGERFS); the halving stop usually ends
 * the steps sooner. More can only help while each step still halves the
 * residual, at a solve and a multiply each. */
#define G_TOL_EPS (0.25)
#define G_REFINE_MAX_STEPS (5)
/* G_FAIL_REL_RES: the solve for g has failed when it returns an error, g is not
 * finite (or the size of its terms overflows), or the residual of the refined g
 * is above this fraction of the size of its terms (as for tol_g): not even its
 * leading digit is right. At the setup of a solve SCS then stops with
 * SCS_FAILED and a message asking for a bug report: that is what catches a
 * broken linear solver. At a scale update SCS instead goes back to the previous
 * R, where g was solved, and continues (scale_update_rollback in src/scs.c; it
 * stops as at setup only if g cannot be solved there either): there the check
 * can also fire on a working solver at the limit of what it can do with the new
 * R, as CG did on beaconfd (10 n iterations without a correct digit, 0.73 after
 * refinement, with R_y spanning 1e9; 3 of 564 starts). The fraction is at most
 * 1 for every g (g = 0 has 1); a solver solving with an R 1e12 times too large
 * gives 1.00 (test_g_accuracy). The refined g of working solvers is far below
 * it, at about a machine epsilon or less (with G_TOL_EPS 4 it was at most 1e-15
 * in double precision and 5e-7 in single, on the sets below) except where the
 * solver itself cannot do better: at most 5e-3 with CG at its iteration limit
 * (QSCFXM1, double precision), 2e-3 with QDLDL in single precision
 * (neos-4413714-turia), and 1.5e-2 with test_refine_perturbed's solver on R
 * (1 - 0.55); first solves, before refinement, reached 0.17 (single-precision
 * CG, neos-3656078-kumeu) and were refined to 5e-7 with G_TOL_EPS 4 (measured
 * with QDLDL and CG in both precisions on the test suite, 250 random SOCPs and
 * 17 badly scaled LPs and QPs, 14 of them with CG in double precision; dense
 * and Accelerate on the test suite). Lower would also stop solvers that are
 * inaccurate rather than broken, which SCS's iterations cope with (and roll
 * back more updates); higher would let through a g with no correct digit. */
#define G_FAIL_REL_RES (0.5)

/* --- Accuracy of each iteration's solve ---
 * Every iteration solves K u = [R_x v_x; -R_y v_y] (K as above), and SCS's
 * residuals cannot fall much below the error of that solve. On iterations
 * that are multiples of REFINE_CHECK_INTERVAL SCS computes the residual of
 * the solve, whatever the linear solver, and compares its inf-norm with tol,
 * the smaller inf-norm of SCS's primal and dual residual vectors at the last
 * convergence check. When the residual is above an "on" factor times tol,
 * refinement turns on; once on, it turns off at a later check whose residual
 * is at or below a smaller "off" factor times tol (with one factor for both
 * it switched off and on again, and the unrefined solves in between undid
 * the progress: at 1e-9 that stalled 20 of 400 runs with cuDSS). In double
 * precision it is never on while the residual is at rounding level, where
 * refinement cannot lower it (REFINE_ROUNDING_EPS). What "on" means depends
 * on the kind of linear solver, chosen at compile time by INDIRECT, so each
 * solver behaves the same in every build and interface:
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
 * the inexact solves themselves. SCS keeps the previous iteration's solution
 * (refine_prev) while refinement is on and passes it as the warm start of
 * scs_solve_lin_sys (linsys.h), which makes the correction form, its
 * fallback and their multiplies by K; SCS decides only when, by the rules
 * here. A linear solver can compute the multiplies with its own copy of K
 * (MKL Pardiso; cuDSS from CUDSS_MIN_PRODUCT_NNZ, where the residual stays
 * on the GPU and the warm start is copied there), else they are SCS's
 * products with the problem data (linsys/scs_matrix.c).
 *
 * Indirect solvers (conjugate gradients, CPU and GPU). CG already solves to a
 * tolerance tied to SCS's residuals, and its error is set by that tolerance,
 * not by the size of the right-hand side, so the correction form would not
 * make it more accurate. What goes wrong is that a solve misses its
 * tolerance: it stops at its iteration limit, the residual it tracks drifts
 * from the true one (in single precision), or y, which CG computes from x,
 * carries rounding that CG's reduced system does not see. While on, each
 * iteration's solve is one call (scs_solve_lin_sys with max_refine_steps
 * REFINE_CG_MAX_STEPS, linsys.h) that computes the true residual
 * r = rhs - K u in the full space (one multiply) and, while ||r||_inf is
 * above the CG tolerance, adds a cold CG solve of K d = r to that tolerance,
 * keeping each step only if it lowers ||r|| and stopping after one that
 * does not halve it. So an iteration with refinement on costs one multiply,
 * and CG iterations only where its solve missed. At a check while on, a
 * plain CG solve of the same right-hand side from the same start, made only
 * for the decision, decides whether to stay on, as for the direct solvers
 * (one more CG solve per REFINE_CHECK_INTERVAL iterations while on), so
 * once on, refinement stays on for at least REFINE_CHECK_INTERVAL
 * iterations. The earlier rule decided at every iteration while on, with
 * SCS's own residual, and added a cold CG solve of r to CG_BEST_TOL each
 * time, whether or not the solve had missed. Against it (KKT-verified at 10
 * times eps; CPU CG with clang on Apple silicon and GCC on ARM; random SOCPs
 * at 1e-9 and, in single precision, 1e-6, the Netlib and SDPLIB certificate
 * problems at 1e-4 and 1e-6, badly scaled QPs and LPs at 1e-4; from one
 * start, and from 8 more for every problem that moved and, with GCC, for
 * all 100 single-precision SOCPs) the verified solves were the same within
 * noise. In double precision 2 of 1,138 runs were lost and none gained:
 * refinery (an infeasible LP) at one of its 9 starts, at both 1e-4 and 1e-6
 * (at 1e-4 its certificate failed the check, at 1e-6 the 120 s limit stopped
 * it). In single precision 48 runs were gained and 48 lost of 2,883 (outcomes
 * there flip both ways with rounding; s30_004 was solved from 5 of its 9
 * starts instead of 8, s30_009 from 8 instead of 6). The time, in shifted
 * geometric mean averaged over the starts, was 0.994 and 0.987 times in
 * single precision on the two machines (1.021 and 1.046 from one start), and
 * 1.004 and 0.999 in double precision, where most problems are random SOCPs
 * on which CG refinement never turns on; on Apple silicon the total time of
 * the runs in double precision, dominated by certificates and problems that
 * stall, was 1.266 times (from one start). Over all 1,792 problems the time
 * was 0.998 times averaged over the starts (1.013 from one start). These
 * figures predate the solve for g reaching machine precision; since then,
 * from one start on 100 of the single-precision SOCPs (Apple silicon), the
 * verified solves were the same (2 gained, none lost) and the time 1.124
 * times (1.068 before), most of it on the 20 largest, which both rules now
 * run to the iteration limit with refinement on for 95% of the iterations.
 * Where CG refinement fires in double precision, on certificates and on
 * problems that stall, it was on for 2.5 to 6 times as many iterations, 40 to
 * 45% of the steps were rejected (no lower residual), and some runs took much
 * longer: QSCFXM2 for a fixed 2,000 iterations 2.3 times the time, refinery
 * at 1e-6 1.9 times (from one start). Deciding at every iteration while on
 * instead (the plain solve with res, then, while on, a second call that
 * refines that solution) kept double precision as it was (on for as many
 * iterations as the earlier rule, no run lost; QSCFXM2 1.6 times, refinery
 * 1.5, the total time on Apple silicon 1.058 times), but in single precision,
 * where refinement is on for most of the solve, it took 2.5 to 4% more CG
 * iterations and more time, with the same solves within noise: 1.053 and
 * 1.018 times averaged over the starts on the two machines, against 0.994 and
 * 0.987 (over all problems 1.026 against 0.998), so refinement is decided at
 * checks. With GPU CG (a GB10) in single precision, on 20 SOCPs from one
 * start, 20 were solved against 19, at 1.06 times the time; in double
 * precision, on 10 of them at 1e-9, refinement never turned on. A step must
 * solve for the residual in the full space: CG restarted on rhs from u, the
 * same in exact arithmetic, works on CG's reduced residual and recomputes y
 * from x, so it cannot remove the rounding the full residual measures (in
 * single precision, from the same u, it lowered the residual by 10% or more
 * in 38% of 455,177 steps, against 96% for the cold step, and of 25 random
 * SOCPs at 1e-5 it solved 16 instead of 25, in 12 times the time).
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
 * of 78 runs). In single precision it tends to stay on once it turns on
 * (REFINE_ROUNDING_EPS). */
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
 * iteration costs one more multiply, and CG solves where it missed its
 * tolerance. A smaller factor refines more (measured with the earlier rule,
 * a CG solve every iteration while on): 0.1 took slightly fewer
 * iterations at 1e-9 with QDLDL but fired in 15 of the 336 runs at 1e-4 and
 * 1e-6 instead of 9, and where refinement fires it cost more (3% more
 * linear-system time with CG). A larger one leaves solves that limit
 * accuracy unrefined. */
#define REFINE_INDIRECT_ON_FACTOR (1.0)
/* REFINE_INDIRECT_OFF_FACTOR: once on, CG refinement stays on until the
 * solve's error is at or below this fraction of tol. With an off factor of 1
 * the residuals fell to the error of the refined solves, refinement turned
 * off, and until the next check the unrefined solves undid the progress.
 * The price is that refinement stays on longer once it fires: with CG it
 * refined about twice as many solves (this paragraph was measured with the
 * earlier rule, a CG solve every iteration while on). The off threshold is
 * at least CG_BEST_TOL, the smallest tolerance CG is given (refine_next_on
 * in src/scs.c: tol is at least CG_BEST_TOL / REFINE_INDIRECT_OFF_FACTOR,
 * 1e-5 in single precision, 1e-11 in double), which in single precision is
 * near what a refined solve can reach, so this happens more often there: at
 * the default tolerances, on 30 LPs and QPs from 3 starts each, refinement
 * was on for most of the solve in 41 of 90 runs with CG (21 with an off
 * factor of 1), mostly runs that stall, where it cost up to 2.1 times the
 * linear-system time per iteration. It also solved more: 47 instead of 41
 * (gaining 8, losing 2). */
#define REFINE_INDIRECT_OFF_FACTOR (0.1)
/* REFINE_ROUNDING_EPS: in double precision the per-iteration refinement
 * turns on, and stays on, only while the plain solve's residual
 * ||rhs - K u||_inf is above REFINE_ROUNDING_EPS machine epsilons
 * (SCS_MACH_EPS) times || |rhs| + |K| |u| ||_inf, besides the factors times
 * tol above (refine_next_on in src/scs.c). That size is, row by row, the sum
 * of the magnitudes of the terms the residual is computed from, so the
 * rounding of the computed residual is a few epsilons of it whatever the
 * scaling of K, rhs and u (the error bound of a sum of products): below it
 * refinement can lower ||r||_inf by at most a small factor (to about one
 * epsilon of that size, the rounding of the row with the largest terms),
 * and the computed residual cannot show even that. It replaces an absolute
 * floor on tol of 1e-12 (an earlier rule), which was below rounding level
 * where those terms are large: where tol goes to zero (one of SCS's
 * residuals on a certificate) refinement turned on and stayed on to the end,
 * chasing rounding (test_refine_certificate, 1000 iterations: on for 799
 * with QDLDL on an unbounded LP and 703 with CG on an infeasible one; now
 * 0). 4: with QDLDL (random SOCPs at 1e-9) the correction form's residuals
 * were at most 1.6 epsilons of the size of their terms in 95% of the checks
 * (median 0.5; the rest had not converged), and plain solves on
 * certificates 0.6 to 2.7, so a plain solve within 4 is as accurate as
 * refinement could make it; 2 left refinement on, or switching, for up to
 * 300 of 500 iterations on certificates; a larger value turns off refinement
 * that still lowers the residual severalfold. On 500 random SOCPs at 1e-4,
 * 1e-6 and 1e-9 (QDLDL, dense, Accelerate), on the Netlib and SDPLIB
 * certificate problems at 1e-4 and 1e-6 and on badly scaled QPs and LPs at
 * 1e-4 and 1e-6 (regress340; QDLDL, KKT-verified at 10 times eps) every run
 * was solved as with the absolute floor, and every run took the same
 * iterations but two (and two stopped by a time limit in both): LASER at
 * 1e-4, where tol fell below 1e-14 while the plain residual, 1e-10 of its
 * terms, was far above their rounding, so refinement now turned on (for 11%
 * of the iterations; 0.7% fewer iterations), and klein1 (an infeasible LP)
 * at 1e-6, where the plain residual was 1.7 epsilons of its terms and
 * refinement stayed off where the floor let it on (2,100 iterations instead
 * of 1,150; from 8 more starts, equal from 6 and more from 2). The suite's
 * certificate problems mostly find their certificate before tol reaches
 * rounding level, so this rule matters where SCS keeps iterating on one.
 * The size is row by row: ||rhs|| alone misses the cancellation in K u (the
 * row-by-row size was up to 150 times it), and the normwise ||rhs|| + ||K||
 * ||u|| was up to 1e3 times the row-by-row size on the random problems and a
 * median 1e5 (up to 2.5e6) on regress340, where it would stop refinement that
 * is needed. Computing it costs a pass over A and P, made only at a check
 * where it can change the outcome. Single precision has no such test (0), and
 * keeps refinement on as before: there the correction form, and CG's
 * refinement step, also lower an error of the solution that the residual does
 * not show (CG and the dense solver's reduced systems: plain residuals of 0.1
 * to 0.5 epsilons of the size of their terms while refinement moved the
 * solution by 6 to 17 epsilons of ||u||, and runs stalled without it). With
 * the test at 4 in single precision, at 1e-6 (KKT-verified at 1e-5): CG
 * solved 94 and 97 of 100 random SOCPs on two machines, against 100 and 98
 * without it (the earlier CG rule, see "Indirect solvers" above: 98 and 99),
 * and the dense solver 489 of 500 against 497 (with the test on staying on
 * only: 496 against 497, and Accelerate 493 against 494); at 1e-5 nothing
 * changed but the time. The price is that in single precision refinement
 * stays on once the plain solves are at rounding level, as on certificates.
 * The absolute floor's other role stays: with CG tol is at least CG_BEST_TOL
 * / REFINE_INDIRECT_OFF_FACTOR, the old floor in single precision (1e-5: the
 * single-precision solves that need refinement have tol near 1e-3, and 1e-6
 * refined more without converging faster, test_adaptive_scale with CG and an
 * off factor of 1: 54 refinement steps and 225 iterations, against 22 and
 * 175) and 1e-11 in double. The direct solvers need no floor in single
 * precision: without the old one every run of 500 random SOCPs at 1e-6 took
 * the same iterations with QDLDL, dense and Accelerate. */
#ifdef SFLOAT
#define REFINE_ROUNDING_EPS (0.)
#else
#define REFINE_ROUNDING_EPS (4.)
#endif
/* REFINE_CG_MAX_STEPS: while CG refinement is on, the most refinement steps
 * an iteration's call takes (scs_solve_lin_sys's max_refine_steps). Each is
 * a cold CG solve, to the iteration's CG tolerance, for the residual in the
 * full space, and a multiply by K, taken only while that residual is above
 * the tolerance; a step is kept only if it lowers the residual, and the
 * steps stop after one that does not halve it (LIN_SYS_REFINE_CONTRACTION).
 * In single precision at 1e-6, where refinement fires most (random SOCPs),
 * 99.6% of the calls took a step (CG's tolerance, at its floor CG_BEST_TOL
 * there, is below the rounding of the residual), 1.5 to 1.6 steps per call,
 * and 83 to 84% of the steps were kept (since the solve for g reaches
 * machine precision, on 100 of them from one start: 99.9%, 1.8 steps and
 * 81%). One step solved 97 of 100 SOCPs
 * where two solved 100 (the earlier CG rule: 98), at 0.94 times the time of
 * two (clang on Apple silicon); more steps add CG solves that the halving
 * stop mostly ends. In double precision, where CG refinement fires on
 * certificates and on problems that stall, calls took about 1.5 steps and
 * 40 to 45% of the steps were rejected (no lower residual), each costing a
 * CG solve. One step in double precision (two in single) solved the same
 * within noise (one run fewer of 1,138: refinery at 1e-6 from one start,
 * stopped by the 120 s limit), at 1.006 and 1.015 times the earlier rule's
 * time in shifted geometric mean averaged over the starts (two steps: 1.004
 * and 0.999), though QSCFXM2 for a fixed 2,000 iterations took 1.6 times the
 * earlier rule's time instead of 2.3. */
#define REFINE_CG_MAX_STEPS (2)

/* CUDSS_MIN_PRODUCT_NNZ: the smallest K (nonzeros of the stored triangle,
 * diagonal included) for which the cuDSS solver computes the refinement's
 * multiplies by K on the GPU, with its own copy of K; below it they are SCS's
 * products on the CPU, from the problem data (linsys/scs_matrix.c), as for the
 * solvers without products of their own. It is here rather than in the cuDSS
 * solver so that the tests see the same value (they check which K a residual is
 * measured with). Set it with make CUDSS_MIN_PRODUCT_NNZ=... (scs.mk passes it
 * to the C and CUDA compilers alike), or with CMake
 * -DCMAKE_C_FLAGS=-DCUDSS_MIN_PRODUCT_NNZ=...; 0 computes them on the GPU at
 * every size (a CI job does, so that path is tested). A GPU product costs about
 * 0.04 ms whatever the size (the copy of b, kernel launches, the amax and its
 * scalar copy), SCS's own about 1 to 1.5 ns per nonzero. Measured per product
 * on an NVIDIA GB10 (Grace CPU, unified memory, 33 problems): the GPU was 1.3
 * to 50 times slower up to 4.2e4 nonzeros (0.04 against 0.0008 ms at 1.2e3) and
 * faster from 9.0e4 on (0.64 to 0.87 times up to 1.2e5, 0.12 to 0.37 times from
 * 3.7e5; BOYD1, 9.3e5: 0.18 against 1.27 ms). In the solves, on random LPs and
 * SOCPs with 1.3e4 nonzeros at 1e-9, an iteration with the correction form on
 * cost 1.6 times a plain one with the GPU products against 1.1 with SCS's
 * (measured while the warm start stayed on the GPU; it is now copied there from
 * the host, which on 10 larger problems added 2 to 9% of a plain solve to such
 * an iteration). A discrete GPU, where copies cross PCIe, likely moves the
 * break-even point up; not measured. Deterministic: the choice depends on the
 * problem's size only. At or above it the products are part of the solver: if
 * they cannot be set up (GPU memory for the second copy of K, a cuBLAS or
 * cuSPARSE error, K too large for the index type), scs_init_lin_sys_work fails
 * with a message, rather than leaving the products to SCS. Setup cost (NVIDIA
 * GB10, random K): 3, 35 and 340 ms at 1e5, 1e6 and 8e6 nonzeros, 5 to 13% of
 * cuDSS's analysis and factorization; of that cusparseCreate and cublasCreate
 * took 0.4 to 0.8 ms. */
#ifndef CUDSS_MIN_PRODUCT_NNZ
#define CUDSS_MIN_PRODUCT_NNZ (60000)
#endif

/* MKL_NO_OWN_PRODUCTS: 1 makes the MKL Pardiso solver compute the
 * refinement's multiplies by K with SCS's products (linsys/scs_matrix.c) at
 * every size instead of MKL's sparse BLAS, so that this path, otherwise taken
 * only for a K beyond MKL's 32-bit sparse integers, stays tested (a CI job
 * builds with it). It is slower: an iteration with the correction form on
 * costs 1.45 instead of 1.32 times a plain one (33 LPs and QPs at 1e-9, 4
 * threads, AMD Zen 3). 0, the default, uses MKL's. It is here so that the
 * tests see the same value. Set it with make MKL_NO_OWN_PRODUCTS=1 or CMake
 * -DCMAKE_C_FLAGS=-DMKL_NO_OWN_PRODUCTS=1; it is on only for a value above 0,
 * as NO_PRINTING and the other switches tested with #if. */
#ifndef MKL_NO_OWN_PRODUCTS
#define MKL_NO_OWN_PRODUCTS (0)
#endif

/* --- Conjugate gradient (CG) parameters, only used with indirect solver --- */
/* Accuracy of the cold CG solve for g0 (see "Accuracy of the solve for g"
 * above; its refinement then goes to about machine precision) and the floor
 * of the per-iteration (warm) ones, set in src/scs.c.
 * Single precision cannot reach 1e-12: a CG asked for it
 * runs to its iteration limit on every solve and its rounding errors then
 * reach the iterates (issue_499's infeasibility certificate stalled with
 * tau at rounding level). 1e-6 is about ten times float's machine epsilon;
 * 1e-7 also passes the test suite, 1e-5 starts losing accuracy. */
#ifdef SFLOAT
#define CG_BEST_TOL (1e-6)
#else
#define CG_BEST_TOL (1e-12)
#endif
/* CG_ZERO_RHS_TOL: a CG solve returns x = 0 without iterating when
 * ||b||_inf is at most MIN(tol, CG_ZERO_RHS_TOL); x = 0 then meets tol, its
 * residual being ||b||. The per-iteration solves have tol >= CG_BEST_TOL, so
 * for them this is the 1e-12 of SCS 3.3; the MIN honours a smaller tol (the
 * solve for g and its refinement steps). Higher returns 0 for right-hand
 * sides that a solve would still resolve; lower runs CG on right-hand sides
 * at rounding level. */
#define CG_ZERO_RHS_TOL (1e-12)
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
/* eigCG deflation in the CPU indirect solver (Stathopoulos & Orginos
 * 2010). The two settings below are fixed when SCS is compiled, so a
 * build has one behaviour and nothing is read at run time. Override them
 * with -D: make CFLAGS="-DDEFLATE_VECTORS=0" with the Makefile, or
 * -DCMAKE_C_FLAGS="-DDEFLATE_VECTORS=0" when configuring with CMake.
 * Deflation needs USE_LAPACK and double precision: without LAPACK or with
 * SFLOAT it is off whatever the values. The GPU indirect solver does not
 * deflate. */

/* DEFLATE_VECTORS: number of approximate eigenvectors, for the smallest
 * eigenvalues of the preconditioned reduced matrix, harvested from the CG
 * iterations of each cold solve (g = K^{-1}h, once per metric change) and
 * used to deflate the warm solves that follow. Cost: 2 * n floats per
 * vector kept in memory, one dot product and two axpys of length n per
 * vector in each solve (no matrix multiplies), and one CG iteration's
 * matrix multiplies per vector at each harvest.
 * Tradeoff: more vectors take more of the small end of the spectrum out
 * of the warm solves, which then need fewer CG iterations on
 * ill-conditioned problems, for more memory and more work per solve and
 * per harvest. 0 turns deflation off (plain preconditioned CG, no extra
 * memory). At run time it is capped at n. */
#ifndef DEFLATE_VECTORS
#define DEFLATE_VECTORS (30)
#endif
#if DEFLATE_VECTORS < 0
#error "DEFLATE_VECTORS must be >= 0 (0 turns eigCG deflation off)"
#endif

/* EIGCG_WINDOW: number of CG (Lanczos) vectors a harvest keeps, n floats
 * each, allocated for the cold solve and freed after it. When the window
 * is full, eigCG compresses it to at most 2 * DEFLATE_VECTORS Ritz vectors
 * (a thick restart: dense eigendecompositions of the window's size and an
 * n x window by window x 2 * DEFLATE_VECTORS product). Tradeoff: a larger
 * window restarts less often and keeps more of the CG history between
 * restarts, for more memory during the harvest and larger restarts. At
 * run time the window is raised to at least 3 * DEFLATE_VECTORS and
 * 2 * DEFLATE_VECTORS + 2, so each cycle has room for DEFLATE_VECTORS new
 * vectors: with the default 30 vectors the window is 90, so a value below
 * 90 has no effect unless DEFLATE_VECTORS is lowered as well (the default
 * 60 takes effect only when fewer than 20 vectors are used). Unused when
 * DEFLATE_VECTORS is 0. */
#ifndef EIGCG_WINDOW
#define EIGCG_WINDOW (60)
#endif
#if EIGCG_WINDOW < 0
#error "EIGCG_WINDOW must be >= 0"
#endif

#ifdef __cplusplus
}
#endif
#endif
