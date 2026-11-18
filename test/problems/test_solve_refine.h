#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <string.h>

/*
 * Refinement of each iteration's solve (solve_direct and refine_iter in
 * src/scs.c; "Accuracy of each iteration's solve" in glbopts.h), on an LP
 * with a second-order cone built like the MATLAB interface's mixed-cone test,
 * at a larger size: A dense, 130 x 100, entries uniform in [-1, 1];
 * b = A x + s and c = -A' s for a random x and an s inside the cones (50 LP
 * rows, then one second-order cone of size 80), so that both the primal and
 * the dual are strictly feasible.
 *
 * test_solve_refine_tight solves it at a tolerance where the error of a QDLDL
 * solve is as large as the residuals SCS has to reach: 1e-9 in double
 * precision, the default tolerances in single precision. Without
 * per-iteration refinement QDLDL does not converge there: it reaches
 * max_iters with residuals near 1e-4 in double precision (with and without
 * FMA contraction) and near 0.1 in single precision, where 49 of 50 seeds
 * fail. With it, QDLDL refines (in correction form, the direct solvers' way)
 * and solves it: in 175 iterations in double precision (225 with one
 * refinement step per iteration instead, the rule before the correction
 * form). CG, the dense solver and Apple Accelerate solve this system
 * accurately enough that refinement need not fire there, so only QDLDL is
 * required to refine.
 *
 * test_solve_refine_idle solves it at the default tolerances in double
 * precision: the residual is checked and refinement does not fire, whatever
 * the linear solver.
 *
 * test_solve_refine_csv repeats the tight solve with QDLDL in double precision
 * while logging to CSV, which recomputes the residuals every iteration, and
 * requires the same refinement (iterations, residual checks and refined
 * iterations) and the same objectives to 1e-6.
 *
 * "refined" below counts the refinement steps (indirect solvers) plus the
 * iterations solved in correction form (direct solvers); a build has only
 * one of the two.
 */

/* the instance (seeds 10 and 43 also pass in double precision; in single
 * precision seed 43 sometimes converges without refining, in 7 of 404 runs
 * from perturbed initial scales) */
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
 * checks and refined solves made (refinement steps plus iterations solved in
 * correction form) */
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
    *steps = w->refine_steps + w->refine_corrections;
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
             "checks, %li refined\n",
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
             "checks, %li refined\n",
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
 * check only, so logging must not change when it fires: the tight solve with
 * and without logging must refine the same number of solves. Comparing with
 * residuals recomputed at the solve instead doubled them (50 iterations in
 * correction form instead of 25, in 175 iterations either way; with the rule
 * before the correction form it also changed the iterations, 200 instead of
 * 225). Only that count is compared exactly: from 101 initial scales within
 * 50 ulps of the default it is 25 in every run, with and without FMA
 * contraction, while the iterations (175 or 200) and the residual checks (7
 * or 8) move with the last bits, and with USE_LAPACK every vector operation
 * goes through BLAS, where MKL's BLAS can round differently when the logged
 * solve's arrays are aligned differently in memory (it did in CI). The
 * objectives are compared to 1e-6 for the same reason. In single precision
 * the count is not robust either (the rule before the correction form gave
 * 26 refinement steps from 100 of 101 scales and 1 from the other), so the
 * test is skipped there; the bug it guards against does not depend on the
 * precision. Checked with QDLDL, where refinement fires on this instance
 * (the CG solvers' own tolerance follows the logged residuals, so there
 * logging changes the iterates whether or not refinement fires). */
static const char *test_solve_refine_csv(void) {
#ifdef SFLOAT
  return mu_skipped;
#else
  scs_int status0, status1, checks0, checks1, steps0, steps1;
  ScsInfo info0, info1;
  if (!strstr(scs_get_lin_sys_method(), "qdldl")) {
    return mu_skipped;
  }
  status0 = solve_refine_solve(SOLVE_REFINE_TIGHT_EPS, SCS_NULL, &info0,
                               &checks0, &steps0);
  status1 = solve_refine_solve(SOLVE_REFINE_TIGHT_EPS, "test_solve_refine.csv",
                               &info1, &checks1, &steps1);
  remove("test_solve_refine.csv");
  scs_printf("test_solve_refine_csv: without logging %li iterations, %li "
             "checks, %li refined, pobj %.17g; with logging %li iterations, "
             "%li checks, %li refined, pobj %.17g\n",
             (long)info0.iter, (long)checks0, (long)steps0, info0.pobj,
             (long)info1.iter, (long)checks1, (long)steps1, info1.pobj);
  mu_assert("test_solve_refine_csv: not solved",
            status0 == SCS_SOLVED && status1 == SCS_SOLVED);
  mu_assert("test_solve_refine_csv: logging changed the refinement",
            steps0 == steps1);
  mu_assert("test_solve_refine_csv: logging changed the result",
            ABS(info0.pobj - info1.pobj) <= 1e-6 * MAX(1., ABS(info0.pobj)) &&
                ABS(info0.dobj - info1.dobj) <=
                    1e-6 * MAX(1., ABS(info0.dobj)));
  return SCS_NULL;
#endif
}

/*
 * test_solve_refine_correction drives the direct solvers' correction form on
 * every direct solver and checks the solution against the problem data. The
 * instance is the one above with its second column replaced by the first
 * plus SOLVE_REFINE_DEP times the second (and b, c rebuilt so the primal and
 * dual stay strictly feasible): the near-dependent columns make the KKT
 * matrix ill-conditioned enough that every direct solver turns refinement on
 * at 1e-9 (QDLDL, the dense solver and Accelerate checked from 21 initial
 * scales within 10 ulps of the default: every run solved, solving 74 to 749
 * of its iterations in correction form; MKL Pardiso and cuDSS as well). In
 * single precision it runs at the default tolerances.
 *
 * Each direct solver solves it as set up and, if the linear solver computes
 * the refinement's products, again with SCS computing them itself
 * (refine_hook = 0). Every run must be solved with a solution whose relative
 * primal and dual residuals and gap, recomputed from the data, are at most
 * 1e-5 (1e-2 in single precision): a correction-form bug gives a wrong
 * answer of order one or no convergence. Only QDLDL, which SCS bundles, is
 * also required to have solved iterations in correction form: whether a
 * third-party solver (dense LAPACK, Accelerate, MKL Pardiso, cuDSS) turns
 * refinement on depends on the accuracy of its solves, which can change
 * with its version, and on iteration counts that change chaotically with
 * rounding; for them the count is printed. The objectives of the two runs
 * are only printed: the products round differently, and how a run continues
 * from there is not reproducible.
 * The indirect solvers skip it: their rule is separate, and without
 * deflation the GPU conjugate gradient solver needs minutes on this
 * ill-conditioned instance at 1e-9.
 */
#define SOLVE_REFINE_DEP (1e-5)
#ifdef SFLOAT
#define SOLVE_REFINE_SOL_TOL (1e-2)
#else
#define SOLVE_REFINE_SOL_TOL (1e-5)
#endif

static void solve_refine_dep_data(ScsData *d, ScsCone *k) {
  scs_int i, j, n, m, ml = 50;
  scs_float *s;
  solve_refine_data(d, k);
  n = d->n;
  m = d->m;
  /* the cone point s of solve_refine_data: 1 in the LP rows, (t, 0.5, ...)
   * with t twice the norm of the rest in the cone */
  s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  for (i = 0; i < ml; ++i) {
    s[i] = 1.;
  }
  s[ml] = SQRTF((scs_float)(m - ml - 1));
  for (i = ml + 1; i < m; ++i) {
    s[i] = 0.5;
  }
  for (i = 0; i < m; ++i) {
    d->A->x[m + i] = d->A->x[i] + SOLVE_REFINE_DEP * d->A->x[m + i];
  }
  /* b = A x + s with x = 0.3 * ones, c = -A' s */
  memcpy(d->b, s, m * sizeof(scs_float));
  memset(d->c, 0, n * sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      d->b[i] += d->A->x[j * m + i] * 0.3;
      d->c[j] -= d->A->x[j * m + i] * s[i];
    }
  }
  scs_free(s);
}

/* the largest of the relative primal residual, dual residual and gap of sol,
 * from the data */
static scs_float solve_refine_kkt_error(const ScsData *d,
                                        const ScsSolution *sol) {
  scs_int n = d->n, m = d->m, i;
  scs_float *pr = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *dr = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float ep, ed, eg, cx, by;
  SCS(accum_by_a)(d->A, sol->x, pr);
  SCS(accum_by_atrans)(d->A, sol->y, dr);
  for (i = 0; i < m; ++i) {
    pr[i] += sol->s[i] - d->b[i];
  }
  for (i = 0; i < n; ++i) {
    dr[i] += d->c[i];
  }
  cx = SCS(dot)(d->c, sol->x, n);
  by = SCS(dot)(d->b, sol->y, m);
  ep = SCS(norm_inf)(pr, m) / (1. + SCS(norm_inf)(d->b, m));
  ed = SCS(norm_inf)(dr, n) / (1. + SCS(norm_inf)(d->c, n));
  eg = ABS(cx + by) / (1. + ABS(cx) + ABS(by));
  scs_free(pr);
  scs_free(dr);
  return MAX(ep, MAX(ed, eg));
}

/* solves the instance, with SCS's own products if own_products; returns the
 * status, the iterations solved in correction form and the solution's error
 * (solve_refine_kkt_error) */
static scs_int solve_refine_dep_solve(scs_int own_products, ScsInfo *info,
                                      scs_int *corrections, scs_int *hook,
                                      scs_float *err) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsWork *w;
  scs_int status = SCS_FAILED;
  memset(info, 0, sizeof(ScsInfo));
  *corrections = *hook = -1;
  *err = -1.;
  solve_refine_dep_data(d, k);
  scs_set_default_settings(stgs);
  if (SOLVE_REFINE_TIGHT_EPS > 0) {
    stgs->eps_abs = SOLVE_REFINE_TIGHT_EPS;
    stgs->eps_rel = SOLVE_REFINE_TIGHT_EPS;
  }
  stgs->max_iters = 100000;
  stgs->verbose = 0;
  w = scs_init(d, k, stgs);
  if (w) {
    if (own_products) {
      w->refine_hook = 0;
    }
    *hook = w->refine_hook;
    status = scs_solve(w, sol, info, 0);
    *corrections = w->refine_corrections;
    scs_finish(w);
    *err = solve_refine_kkt_error(d, sol);
  }
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  return status;
}

static const char *test_solve_refine_correction(void) {
  scs_int status, corrections, hook, own, must_refine;
  scs_float err;
  ScsInfo info;
  if (strstr(scs_get_lin_sys_method(), "indirect")) {
    return mu_skipped;
  }
  /* only the bundled QDLDL must refine, see above */
  must_refine = strstr(scs_get_lin_sys_method(), "qdldl") != SCS_NULL;
  for (own = 0, hook = 1; own <= 1 && hook; ++own) {
    status = solve_refine_dep_solve(own, &info, &corrections, &hook, &err);
    scs_printf("test_solve_refine_correction: %s, products by %s, status "
               "%li, %li iterations, %li in correction form, pobj %.10e, "
               "error %.2e\n",
               scs_get_lin_sys_method(), hook ? "the linear solver" : "SCS",
               (long)status, (long)info.iter, (long)corrections, info.pobj,
               (double)err);
    mu_assert("test_solve_refine_correction: not solved",
              status == SCS_SOLVED);
    mu_assert("test_solve_refine_correction: wrong solution",
              err >= 0 && err <= SOLVE_REFINE_SOL_TOL);
    if (must_refine) {
      mu_assert("test_solve_refine_correction: refinement did not turn on",
                corrections > 0);
    }
  }
  return SCS_NULL;
}
