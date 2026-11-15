#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <string.h>

/*
 * Refinement of each iteration's solve (refine_iter in src/scs.c,
 * REFINE_CHECK_INTERVAL and REFINE_FACTOR in glbopts.h), on an LP with a
 * second-order cone built like the MATLAB interface's mixed-cone test, at a
 * larger size: A dense, 130 x 100, entries uniform in [-1, 1]; b = A x + s
 * and c = -A' s for a random x and an s inside the cones (50 LP rows, then
 * one second-order cone of size 80), so that both the primal and the dual
 * are strictly feasible.
 *
 * test_solve_refine_tight solves it at a tolerance where the error of a QDLDL
 * solve is as large as the residuals SCS has to reach: 1e-9 in double
 * precision, the default tolerances in single precision. Without
 * per-iteration refinement QDLDL does not converge there: it reaches
 * max_iters with residuals near 1e-4 in double precision (with and without
 * FMA contraction) and near 0.1 in single precision, where 49 of 50 seeds
 * fail. With it, QDLDL refines and solves it in every build checked (with and
 * without FMA contraction, with and without Anderson acceleration, each from
 * 101 initial scales up to 50 ulps apart): in 225 iterations in double
 * precision, in at most 225 in single. CG, the dense solver and Apple
 * Accelerate solve this system accurately enough that refinement does not
 * fire there, so only QDLDL is required to refine.
 *
 * test_solve_refine_idle solves it at the default tolerances in double
 * precision: the residual is checked and refinement does not fire, whatever
 * the linear solver.
 *
 * test_solve_refine_csv repeats the tight solve with QDLDL while logging to
 * CSV, which recomputes the residuals every iteration, and requires the same
 * refinement and the same result.
 */

/* the instance (seeds 10 and 43 also pass in double precision, but refine
 * less reliably in single precision) */
#define SOLVE_REFINE_SEED (9)

/* the tight tolerance (0: SCS's defaults) */
#ifdef SFLOAT
#define SOLVE_REFINE_TIGHT_EPS (0.)
#else
#define SOLVE_REFINE_TIGHT_EPS (1e-9)
#endif

/* fills d and k with the instance */
static void solve_refine_data(ScsData *d, ScsCone *k) {
  scs_int n = 100, ml = 50, mq = 80, m = ml + mq, i, j;
  scs_float *x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->m = m;
  d->n = n;
  d->A = A;
  d->b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  A->m = m;
  A->n = n;
  A->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  A->i = (scs_int *)scs_calloc(n * m, sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(n * m, sizeof(scs_float));
  ran_start(SOLVE_REFINE_SEED);
  for (j = 0; j < n; ++j) {
    A->p[j + 1] = (j + 1) * m;
    for (i = 0; i < m; ++i) {
      A->i[j * m + i] = i;
      A->x[j * m + i] = rand_scs_float();
    }
  }
  for (j = 0; j < n; ++j) {
    x[j] = rand_scs_float();
  }
  /* s: 1 in the LP rows; in the cone (t, 0.5, ..., 0.5) with t twice the
   * norm of the rest */
  for (i = 0; i < ml; ++i) {
    s[i] = 1.;
  }
  s[ml] = SQRTF((scs_float)(mq - 1));
  for (i = ml + 1; i < m; ++i) {
    s[i] = 0.5;
  }
  /* b = A x + s, c = -A' y with y = s */
  memcpy(d->b, s, m * sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      d->b[i] += A->x[j * m + i] * x[j];
      d->c[j] -= A->x[j * m + i] * s[i];
    }
  }
  k->l = ml;
  k->qsize = 1;
  k->q = (scs_int *)scs_calloc(1, sizeof(scs_int));
  k->q[0] = mq;
  scs_free(x);
  scs_free(s);
}

/* solves the instance at eps (eps <= 0: the default tolerances), logging to
 * csv unless it is SCS_NULL; returns the status, SCS's info, and the residual
 * checks and refinement steps made */
static scs_int solve_refine_solve(scs_float eps, const char *csv,
                                  ScsInfo *info, scs_int *checks,
                                  scs_int *steps) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsWork *w;
  scs_int status = SCS_FAILED;
  memset(info, 0, sizeof(ScsInfo));
  info->iter = *checks = *steps = -1;
  solve_refine_data(d, k);
  scs_set_default_settings(stgs);
  if (eps > 0) {
    stgs->eps_abs = eps;
    stgs->eps_rel = eps;
  }
  stgs->max_iters = 100000;
  stgs->verbose = 0;
  stgs->log_csv_filename = csv;
  w = scs_init(d, k, stgs);
  if (w) {
    status = scs_solve(w, sol, info, 0);
    *checks = w->refine_checks;
    *steps = w->refine_steps;
    scs_finish(w);
  }
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  return status;
}

/* solved with every linear solver, and with QDLDL only because refinement
 * fired */
static const char *test_solve_refine_tight(void) {
  scs_int status, checks, steps;
  ScsInfo info;
  status = solve_refine_solve(SOLVE_REFINE_TIGHT_EPS, SCS_NULL, &info, &checks,
                              &steps);
  scs_printf("test_solve_refine_tight: %s, status %li, %li iterations, %li "
             "checks, %li refinement steps\n",
             scs_get_lin_sys_method(), (long)status, (long)info.iter,
             (long)checks, (long)steps);
  mu_assert("test_solve_refine_tight: not solved", status == SCS_SOLVED);
  if (strstr(scs_get_lin_sys_method(), "qdldl")) {
    mu_assert("test_solve_refine_tight: QDLDL solved without refinement, "
              "which this instance does not allow",
              steps > 0);
  }
  return SCS_NULL;
}

/* At the default tolerances the solve is accurate enough: the residual is
 * checked and refinement does not fire, with every linear solver. In single
 * precision the default tolerances are the tight case above. */
static const char *test_solve_refine_idle(void) {
#ifdef SFLOAT
  return mu_skipped;
#else
  scs_int status, checks, steps;
  ScsInfo info;
  status = solve_refine_solve(0., SCS_NULL, &info, &checks, &steps);
  scs_printf("test_solve_refine_idle: %s, status %li, %li iterations, %li "
             "checks, %li refinement steps\n",
             scs_get_lin_sys_method(), (long)status, (long)info.iter,
             (long)checks, (long)steps);
  mu_assert("test_solve_refine_idle: not solved", status == SCS_SOLVED);
  mu_assert("test_solve_refine_idle: the residual was never checked",
            checks > 0);
  mu_assert("test_solve_refine_idle: refinement fired at default tolerances",
            steps == 0);
  return SCS_NULL;
#endif
}

/* Logging to CSV (log_csv_filename) recomputes SCS's residuals on every
 * iteration. Refinement compares with the residuals of the last convergence
 * check only, so logging must not change when it fires or the iterates: the
 * tight solve with and without logging must agree exactly. Checked with QDLDL,
 * where refinement fires on this instance and the solve is deterministic (the
 * CG solvers' own tolerance follows the logged residuals, so there logging
 * changes the iterates whether or not refinement fires). */
static const char *test_solve_refine_csv(void) {
  scs_int status0, status1, checks0, checks1, steps0, steps1;
  ScsInfo info0, info1;
  if (!strstr(scs_get_lin_sys_method(), "qdldl")) {
    return mu_skipped;
  }
  status0 = solve_refine_solve(SOLVE_REFINE_TIGHT_EPS, SCS_NULL, &info0,
                               &checks0, &steps0);
  status1 = solve_refine_solve(SOLVE_REFINE_TIGHT_EPS, "test_solve_refine.csv",
                               &info1, &checks1, &steps1);
  scs_printf("test_solve_refine_csv: without logging %li iterations, %li "
             "checks, %li steps; with logging %li iterations, %li checks, %li "
             "steps\n",
             (long)info0.iter, (long)checks0, (long)steps0, (long)info1.iter,
             (long)checks1, (long)steps1);
  mu_assert("test_solve_refine_csv: not solved",
            status0 == SCS_SOLVED && status1 == SCS_SOLVED);
  mu_assert("test_solve_refine_csv: logging changed the refinement or the "
            "iterates",
            info0.iter == info1.iter && checks0 == checks1 &&
                steps0 == steps1 && info0.pobj == info1.pobj &&
                info0.dobj == info1.dobj && info0.res_pri == info1.res_pri &&
                info0.res_dual == info1.res_dual);
  return SCS_NULL;
}
