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
 *   every direct solver solves with refinement on (refine_corrections > 0),
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
 *   factorization; the test runs RP_DIVERGE_ITERS iterations. Measured in
 *   double precision: 33 to 49 safeguards in 10,000 iterations with QDLDL,
 *   dense, Accelerate and cuDSS below its threshold. A safeguard that never
 *   fires, or that compares the wrong residuals, gives none.
 * - CG solvers, delta = RP_DELTA_CONTRACT, RP_CG_ITERS iterations: the
 *   refinement turns on (refine_steps > 0), which shows that the residual of
 *   the CG solve reaches SCS's on/off rule. CG's refinement step is a second
 *   CG solve warm-started from the first one's solution, which solves with
 *   K' again, so it cannot remove this mismatch (in a real solve CG's K is
 *   SCS's, the same arrays) and the problem is not solved: only the steps
 *   are checked. Measured from 101 initial scales within 50 ulps of the
 *   default: refinement on from iteration 75 in every run, 425 steps in the
 *   500 iterations in double precision at 1e-6, 375 to 400 in single
 *   precision at the default tolerances.
 *
 * Skipped where the linear solver measures its residuals with its own K,
 * which follows the update, so that it sees no mismatch and nothing to
 * refine (the rule each states up front): MKL Pardiso (unless built with
 * MKL_NO_OWN_PRODUCTS), and cuDSS from CUDSS_MIN_PRODUCT_NNZ nonzeros
 * (glbopts.h; this instance has about 1.3e4, so the CI build with the
 * threshold at 0 skips it).
 */

/* RP_DELTA_CONTRACT, RP_DELTA_DIVERGE: the relative changes of R (see
 * above). RP_DIVERGE_ITERS: the iterations of the run that cannot converge,
 * enough for 33 or more safeguards. RP_CG_ITERS: the iterations of the CG
 * run, enough for refinement to turn on (it does from iteration 75) and
 * short, since that run cannot converge either (about 0.1 s natively, and
 * the suite also runs under valgrind). */
#define RP_DELTA_CONTRACT (1e-3)
#define RP_DELTA_DIVERGE (-0.55)
#define RP_DIVERGE_ITERS (10000)
#define RP_CG_ITERS (500)
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
                        scs_int *corrections, scs_int *steps,
                        scs_int *safeguards, scs_float *err) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs_float *dr = SCS_NULL;
  ScsWork *w;
  scs_int status = SCS_FAILED, i;
  memset(info, 0, sizeof(ScsInfo));
  *corrections = *steps = *safeguards = -1;
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
    /* the solver's R, kept alive until scs_finish (CG keeps the pointer) */
    dr = (scs_float *)scs_calloc(d->n + d->m, sizeof(scs_float));
    for (i = 0; i < d->n + d->m; ++i) {
      dr[i] = w->diag_r[i] * (1. + delta);
    }
    if (scs_update_lin_sys_diag_r(w->p, dr) == 0) {
      status = scs_solve(w, sol, info, 0);
      *corrections = w->refine_corrections;
      *steps = w->refine_steps;
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
  scs_int status, corrections, steps, safeguards, cg, nnz;
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
#ifndef MKL_NO_OWN_PRODUCTS
      strstr(method, "mkl") ||
#endif
      (strstr(method, "cuDSS") && nnz >= CUDSS_MIN_PRODUCT_NNZ)) {
    scs_printf("test_refine_perturbed: %s measures with its own K (K has %li "
               "nonzeros), which follows the update: skipped\n",
               method, (long)nnz);
    return mu_skipped;
  }
  cg = strstr(method, "indirect") != SCS_NULL;
  status = rp_solve(RP_DELTA_CONTRACT, cg ? RP_CG_ITERS : 100000, &info,
                    &corrections, &steps, &safeguards, &err);
  scs_printf("test_refine_perturbed: %s, delta %g: status %li, %li "
             "iterations, %li with refinement on, %li refinement steps, %li "
             "safeguards, error %.2e\n",
             method, (double)RP_DELTA_CONTRACT, (long)status, (long)info.iter,
             (long)corrections, (long)steps, (long)safeguards, (double)err);
  if (cg) {
    /* the mismatch cannot be removed by CG (see above); the safeguard is
     * the direct solvers' */
    mu_assert("test_refine_perturbed: CG's refinement did not turn on",
              steps > 0);
    return SCS_NULL;
  }
  mu_assert("test_refine_perturbed: not solved", status == SCS_SOLVED);
  mu_assert("test_refine_perturbed: wrong solution",
            err >= 0 && err <= SOLVE_REFINE_SOL_TOL);
  mu_assert("test_refine_perturbed: refinement did not fire",
            corrections > 0);
  status = rp_solve(RP_DELTA_DIVERGE, RP_DIVERGE_ITERS, &info, &corrections,
                    &steps, &safeguards, &err);
  scs_printf("test_refine_perturbed: %s, delta %g: status %li, %li "
             "iterations, %li with refinement on, %li safeguards\n",
             method, (double)RP_DELTA_DIVERGE, (long)status, (long)info.iter,
             (long)corrections, (long)safeguards);
  mu_assert("test_refine_perturbed: the safeguard did not fire on a "
            "correction form that does not contract",
            safeguards > 0);
  return SCS_NULL;
}
