#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <string.h>

/*
 * The per-iteration refinement (solve_direct and solve_indirect in
 * src/scs.c), and the warm start of scs_solve_lin_sys underneath it, made to
 * fire on purpose. The instance is test_solve_refine's (solve_refine_data,
 * included before this file). After scs_init the linear solver is updated to
 * R' = R (1 + delta) in an array of the test's (scs_update_lin_sys_diag_r),
 * so it solves with K' while SCS's products keep measuring K: a plain
 * solve's residual is about |delta| ||R x||. adaptive_scale is off, so R does
 * not change again.
 *
 * - Direct solvers, contracting, delta = RP_DELTA_CONTRACT: each
 *   correction-form solve shrinks the error it inherits by about delta, so
 *   every direct solver solves with refinement on (refine_on_iters > 0),
 *   and the problem is solved, with a solution whose relative residuals and
 *   gap from the data are at most SOLVE_REFINE_SOL_TOL. A refinement that
 *   does not fire, or a correction form measured with K' (the solver's),
 *   leaves the error at about delta and the solve stalls. Measured in double
 *   precision at 1e-6: 150 iterations, 100 with refinement on, with QDLDL,
 *   dense and Accelerate.
 * - Direct solvers, not contracting, delta = RP_DELTA_DIVERGE: in rows where
 *   R dominates, the correction form multiplies the error it inherits by
 *   about delta / (1 + delta) = -1.2, so it grows from iteration to
 *   iteration until SCS's safeguard replaces it by the plain solve at a
 *   check (refine_safeguards > 0). The problem cannot be solved with such a
 *   factorization; the test runs RP_DIVERGE_ITERS iterations. Nor can the
 *   solve for g be refined to machine precision: its residual is 0.19 of
 *   the size of its terms before refinement and 1.5e-2 after (QDLDL, dense
 *   and Accelerate, double and single precision), which is not a failed
 *   solve (G_FAIL_REL_RES, glbopts.h). Measured in double precision: 29 to
 *   49 safeguards in 10,000 iterations with QDLDL, dense, Accelerate and
 *   cuDSS below its threshold. A safeguard that never fires, or that
 *   compares the wrong residuals, gives none.
 * - CG solvers, delta = RP_DELTA_CONTRACT: CG solves with K'. While
 *   refinement is on, each iteration's solve is one call with refinement
 *   steps (solve_indirect: max_refine_steps REFINE_CG_MAX_STEPS), and each
 *   step solves K' d = r for the residual r that the shared routine computes
 *   with K (SCS's products), which contracts the error by about delta like
 *   the direct solvers' correction form, so the problem is solved with
 *   refinement on (refine_on_iters > 0), to the same accuracy, and the
 *   decisions are made only at checks (refine_checks at most one per
 *   REFINE_CHECK_INTERVAL iterations). A call that takes no steps, or steps
 *   that solve for anything but the full residual (CG restarted on the
 *   right-hand side from x solves with K' again and cannot remove the
 *   mismatch), leaves the error at about delta and the solve stalls.
 *   Measured: 175 iterations, 100 with refinement on in double precision
 *   (75 in single), 7 checks; with the steps removed (max_refine_steps 0
 *   while on) it ran to RP_CONTRACT_ITERS unsolved, error 8e-4, in both.
 *
 * Skipped where the linear solver measures its residuals with its own K,
 * which follows the update, so that it sees no mismatch and nothing to
 * refine (the rule each states up front): MKL Pardiso (unless built with
 * MKL_NO_OWN_PRODUCTS), and cuDSS from CUDSS_MIN_PRODUCT_NNZ nonzeros
 * (glbopts.h; this instance has about 1.3e4, so the CI build with the
 * threshold at 0 skips it).
 */

/* RP_DELTA_CONTRACT, RP_DELTA_DIVERGE: the relative changes of R (see
 * above). RP_CONTRACT_ITERS: the iteration limit of the runs that converge,
 * which take 100 to 175: about ten times that, so that a refinement that
 * does not contract fails in a fraction of a second rather than after
 * 100,000 iterations (the suite also runs under valgrind). RP_DIVERGE_ITERS:
 * the iterations of the run that cannot converge, enough for 29 or more
 * safeguards. */
#define RP_DELTA_CONTRACT (1e-3)
#define RP_DELTA_DIVERGE (-0.55)
#define RP_CONTRACT_ITERS (2000)
#define RP_DIVERGE_ITERS (10000)
#ifdef SFLOAT
#define RP_EPS (0.) /* the default tolerances */
#else
#define RP_EPS (1e-6)
#endif

/* nonzeros of K's stored triangle for this instance (P is NULL), as
 * SCS(form_kkt) builds it */
static scs_int rp_kkt_nnz(const ScsData *d) {
  return d->n + d->m + d->A->p[d->n];
}

/* solves solve_refine_data's instance with the linear solver on R (1 +
 * delta), at most max_iters iterations; returns the status, the counters and
 * the solution's error (solve_refine_kkt_error) */
static scs_int rp_solve(scs_float delta, scs_int max_iters, ScsInfo *info,
                        scs_int *corrections, scs_int *checks,
                        scs_int *safeguards, scs_float *err) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs_float *dr = SCS_NULL;
  ScsWork *w;
  scs_int status = SCS_FAILED, i;
  memset(info, 0, sizeof(ScsInfo));
  *corrections = *checks = *safeguards = -1;
  *err = -1.;
  solve_refine_data(d, k);
  scs_set_default_settings(stgs);
  if (RP_EPS > 0) {
    stgs->eps_abs = RP_EPS;
    stgs->eps_rel = RP_EPS;
  }
  stgs->max_iters = max_iters;
  stgs->verbose = 0;
  stgs->adaptive_scale = 0;
  w = scs_init(d, k, stgs);
  if (w) {
    /* the solver's R, kept alive until scs_finish (CG keeps the pointer).
     * Another array than w->diag_r, on purpose against linsys.h's rule that
     * a new R is written into the array given at init: SCS's products (the
     * shared routine's, and the CG residual in src/scs.c) keep measuring
     * with that array (K) while the solver solves with this one (K'). */
    dr = (scs_float *)scs_calloc(d->n + d->m, sizeof(scs_float));
    for (i = 0; i < d->n + d->m; ++i) {
      dr[i] = w->diag_r[i] * (1. + delta);
    }
    if (scs_update_lin_sys_diag_r(w->p, dr) == 0) {
      status = scs_solve(w, sol, info, 0);
      *corrections = w->refine_on_iters;
      *checks = w->refine_checks;
      *safeguards = w->refine_safeguards;
      *err = solve_refine_kkt_error(d, sol);
    }
    scs_finish(w);
  }
  scs_free(dr);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  return status;
}

static const char *test_refine_perturbed(void) {
  scs_int status, corrections, checks, safeguards, cg, nnz;
  scs_float err;
  ScsInfo info;
  const char *method = scs_get_lin_sys_method();
  {
    ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
    solve_refine_data(d, k);
    nnz = rp_kkt_nnz(d);
    SCS(free_data)(d);
    SCS(free_cone)(k);
  }
  if (
#if MKL_NO_OWN_PRODUCTS <= 0
      strstr(method, "mkl") ||
#endif
      (strstr(method, "cuDSS") && nnz >= CUDSS_MIN_PRODUCT_NNZ)) {
    scs_printf("test_refine_perturbed: %s measures with its own K (K has %li "
               "nonzeros), which follows the update: skipped\n",
               method, (long)nnz);
    return mu_skipped;
  }
  cg = strstr(method, "indirect") != SCS_NULL;
  status = rp_solve(RP_DELTA_CONTRACT, RP_CONTRACT_ITERS, &info, &corrections,
                    &checks, &safeguards, &err);
  scs_printf("test_refine_perturbed: %s, delta %g: status %li, %li "
             "iterations, %li with refinement on, %li checks, %li "
             "safeguards, error %.2e\n",
             method, (double)RP_DELTA_CONTRACT, (long)status, (long)info.iter,
             (long)corrections, (long)checks, (long)safeguards, (double)err);
  mu_assert("test_refine_perturbed: not solved", status == SCS_SOLVED);
  mu_assert("test_refine_perturbed: wrong solution",
            err >= 0 && err <= SOLVE_REFINE_SOL_TOL);
  mu_assert("test_refine_perturbed: refinement did not fire",
            corrections > 0);
  /* Decisions are made only at iterations that are multiples of
   * REFINE_CHECK_INTERVAL, once the first convergence check has set
   * refine_tol (none at iteration 0), so a run of info.iter iterations makes
   * at most info.iter / REFINE_CHECK_INTERVAL. The bound is tight by
   * construction (the CG run: 7 checks in 175 iterations), so a change to
   * when convergence checks run or to how iterations are counted can fail
   * it with no decision between checks; the bound then changes with it. */
  mu_assert("test_refine_perturbed: a decision between checks",
            checks >= 0 && checks <= info.iter / REFINE_CHECK_INTERVAL);
  if (cg) {
    return SCS_NULL; /* the safeguard is the direct solvers' */
  }
  status = rp_solve(RP_DELTA_DIVERGE, RP_DIVERGE_ITERS, &info, &corrections,
                    &checks, &safeguards, &err);
  scs_printf("test_refine_perturbed: %s, delta %g: status %li, %li "
             "iterations, %li with refinement on, %li safeguards\n",
             method, (double)RP_DELTA_DIVERGE, (long)status, (long)info.iter,
             (long)corrections, (long)safeguards);
  mu_assert("test_refine_perturbed: the safeguard did not fire on a "
            "correction form that does not contract",
            safeguards > 0);
  return SCS_NULL;
}
