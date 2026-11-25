#include <float.h>
#include <string.h>

#include "glbopts.h"
#include "minunit.h"
#include "scs.h"
#include "util.h"

/*
 * Shared tiny LP used for all inaccurate termination tests.
 * min x  s.t.  x >= 2
 * n=1, m=1, l=1, A=[[-1]], b=[-2], c=[1], optimal obj=2.
 */
#define _INACCURATE_LP_VARS                                                    \
  scs_float Ax[] = {-1.0};                                                     \
  scs_int Ai[] = {0};                                                          \
  scs_int Ap[] = {0, 1};                                                       \
  scs_float b_feas[] = {-2.0};                                                 \
  scs_float c_feas[] = {1.0}

/*
 * Shared tiny infeasible LP.
 * x >= 1  AND  x <= 0 (infeasible).
 * n=1, m=2, l=2, A=[[-1],[1]], b=[-1,0], c=[1]
 */
#define _INACCURATE_INFEAS_VARS                                                \
  scs_float Ax_i[] = {-1.0, 1.0};                                              \
  scs_int Ai_i[] = {0, 1};                                                     \
  scs_int Ap_i[] = {0, 2};                                                     \
  scs_float b_i[] = {-1.0, 0.0};                                               \
  scs_float c_i[] = {1.0}

/*
 * Shared tiny unbounded LP.
 * min -x  s.t.  x >= 0  (unbounded below).
 * n=1, m=1, l=1, A=[[-1]], b=[0], c=[-1]
 */
#define _INACCURATE_UNBDD_VARS                                                 \
  scs_float Ax_u[] = {-1.0};                                                   \
  scs_int Ai_u[] = {0};                                                        \
  scs_int Ap_u[] = {0, 1};                                                     \
  scs_float b_u[] = {0.0};                                                     \
  scs_float c_u[] = {-1.0}

/*
 * Test SCS_SOLVED_INACCURATE: run a feasible, bounded LP with max_iters=2.
 * After only 2 iterations tau > 0 and kap/tau is small, so set_unfinished
 * returns SCS_SOLVED_INACCURATE.
 */
static const char *test_solved_inaccurate(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  _INACCURATE_LP_VARS;
  scs_int m = 1, n = 1;

  d->m = m; d->n = n; d->b = b_feas; d->c = c_feas;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m; d->A->n = n;
  d->A->x = Ax; d->A->i = Ai; d->A->p = Ap;
  k->l = 1;

  scs_set_default_settings(stgs);
  stgs->max_iters = 2;
  stgs->verbose = 1;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_solved_inaccurate: expected SCS_SOLVED_INACCURATE",
            exitflag == SCS_SOLVED_INACCURATE);
  mu_assert("test_solved_inaccurate: expected max_iters exhausted",
            info.iter == stgs->max_iters);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return 0;
}

/*
 * Test SCS_INFEASIBLE_INACCURATE: run an infeasible LP with eps_infeas=0
 * (disables the infeasibility termination criterion entirely), so the solver
 * runs until max_iters.  set_unfinished then detects kap > tau with the
 * infeasibility direction (bty_tau < 0) and returns SCS_INFEASIBLE_INACCURATE.
 */
static const char *test_infeasible_inaccurate(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  _INACCURATE_INFEAS_VARS;
  scs_int m = 2, n = 1;

  d->m = m; d->n = n; d->b = b_i; d->c = c_i;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m; d->A->n = n;
  d->A->x = Ax_i; d->A->i = Ai_i; d->A->p = Ap_i;
  k->l = 2;

  scs_set_default_settings(stgs);
  stgs->eps_infeas = 0.0; /* disabled: res < 0 never true, forces max_iters */
  /* Use enough iterations so kap >> tau is clearly established across all
   * linear solver backends (numerical trajectory is backend-dependent). */
  stgs->max_iters = 100;
  stgs->verbose = 1;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_infeasible_inaccurate: expected SCS_INFEASIBLE_INACCURATE",
            exitflag == SCS_INFEASIBLE_INACCURATE);
  mu_assert("test_infeasible_inaccurate: expected max_iters exhausted",
            info.iter == stgs->max_iters);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return 0;
}

/*
 * Test SCS_UNBOUNDED_INACCURATE: run an unbounded LP with eps_infeas=0
 * (disables the unboundedness termination criterion entirely), so the solver
 * runs until max_iters.  set_unfinished then detects kap > tau with the
 * unboundedness direction (ctx_tau < 0) and returns SCS_UNBOUNDED_INACCURATE.
 */
static const char *test_unbounded_inaccurate(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  _INACCURATE_UNBDD_VARS;
  scs_int m = 1, n = 1;

  d->m = m; d->n = n; d->b = b_u; d->c = c_u;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m; d->A->n = n;
  d->A->x = Ax_u; d->A->i = Ai_u; d->A->p = Ap_u;
  k->l = 1;

  scs_set_default_settings(stgs);
  stgs->eps_infeas = 0.0; /* disabled: res < 0 never true, forces max_iters */
  /* Use enough iterations so kap >> tau is clearly established across all
   * linear solver backends (numerical trajectory is backend-dependent). */
  stgs->max_iters = 100;
  stgs->verbose = 1;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_unbounded_inaccurate: expected SCS_UNBOUNDED_INACCURATE",
            exitflag == SCS_UNBOUNDED_INACCURATE);
  mu_assert("test_unbounded_inaccurate: expected max_iters exhausted",
            info.iter == stgs->max_iters);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return 0;
}

/*
 * Test SCS_SOLVED_INACCURATE with max_iters=1: the solver exits after a single
 * iteration with tau > 0 (kap/tau small), so set_unfinished returns
 * SCS_SOLVED_INACCURATE.  This exercises the minimum-iteration code path.
 */
static const char *test_max_iters_1(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  _INACCURATE_LP_VARS;
  scs_int m = 1, n = 1;

  d->m = m; d->n = n; d->b = b_feas; d->c = c_feas;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m; d->A->n = n;
  d->A->x = Ax; d->A->i = Ai; d->A->p = Ap;
  k->l = 1;

  scs_set_default_settings(stgs);
  stgs->max_iters = 1;
  stgs->verbose = 1;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_max_iters_1: expected SCS_SOLVED_INACCURATE",
            exitflag == SCS_SOLVED_INACCURATE);
  mu_assert("test_max_iters_1: expected exactly 1 iteration",
            info.iter == 1);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return 0;
}

/*
 * Box LP for the stops at the limit where the last iterate gives no estimate.
 * min 2x  s.t.  -2 <= x <= 2, as x + s_1 = 2, -x + s_2 = 2, s >= 0.
 * n=1, m=2, l=2, A=[[1],[-1]], b=[2,2], c=[2]; solution x = -2, obj -4.
 * SCS solves it in 25 iterations.
 */
static void inaccurate_box_lp(ScsData *d, ScsCone *k) {
  static scs_float Ax[] = {1.0, -1.0};
  static scs_int Ai[] = {0, 1};
  static scs_int Ap[] = {0, 2};
  static scs_float b[] = {2.0, 2.0};
  static scs_float c[] = {2.0};
  d->m = 2;
  d->n = 1;
  d->b = b;
  d->c = c;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = 2;
  d->A->n = 1;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;
  k->l = 2;
}

/* 1 if every entry of x, y and s is NaN */
static int inaccurate_sol_all_nan(const ScsSolution *sol, scs_int n,
                                  scs_int m) {
  scs_int i;
  for (i = 0; i < n; ++i) {
    if (!isnan(sol->x[i])) {
      return 0;
    }
  }
  for (i = 0; i < m; ++i) {
    if (!isnan(sol->y[i]) || !isnan(sol->s[i])) {
      return 0;
    }
  }
  return 1;
}

/*
 * SCS_INDETERMINATE: the box LP stopped at max_iters = 2. After iteration 1
 * the unclipped tau is negative, so project_cones sets tau = 0, and b'y > 0
 * and c'x > 0: the iterate is an estimate of neither a solution (it needs
 * tau > 0) nor a certificate (it needs b'y < 0 or c'x < 0). Measured before
 * this change (QDLDL and CG; double, DLONG and single; normalize,
 * adaptive_scale and acceleration each on or off): unclipped tau -0.023 to
 * -0.046, b'y +0.56 to +0.72 and c'x +0.56 to +1.1 at an iterate of norm
 * about 2, so rounding cannot move the status in either precision. SCS then
 * returned SCS_FAILED, the raw iterate (x = 0.557, which looks like a point)
 * and the status " (inaccurate - reached max_iters)", and it never wrote
 * pobj, dobj, gap, res_pri and res_dual.
 *
 * info starts filled with 0x7f bytes, as from a caller that does not zero it
 * (scs-matlab): before this change the status had no NUL in its 128 bytes,
 * so the suffix was appended past its end (on macOS the fortified strcat
 * aborts the test binary there), and pobj etc. came back as 1.38e306.
 *
 * The other side: max_iters = 3 gives SCS_SOLVED_INACCURATE (tau 0.75, and
 * 0.52 in single precision), the default max_iters gives SCS_SOLVED, and a
 * warm-started solve on the same workspace after SCS_INDETERMINATE solves the
 * problem: the NaN solution is read as zero, a cold start. That solve stops
 * at eps 1e-4 (|x + 2| at most 9e-5 on the Mac builds); a NaN carried into
 * the iterate would give NaN or another status, which the status check and
 * the 1e-2 bound catch.
 */
static const char *test_indeterminate_at_limit(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info;
  ScsWork *w = SCS_NULL;
  scs_int exitflag;
  const char *fail = SCS_NULL;

  inaccurate_box_lp(d, k);
  scs_set_default_settings(stgs);
  stgs->verbose = 1;

  stgs->max_iters = 2;
  memset(&info, 0x7f, sizeof(info));
  exitflag = scs(d, k, stgs, sol, &info);
  if (exitflag != SCS_INDETERMINATE || info.status_val != SCS_INDETERMINATE) {
    fail = "test_indeterminate_at_limit: max_iters 2: expected "
           "SCS_INDETERMINATE";
  } else if (memchr(info.status, '\0', sizeof(info.status)) == SCS_NULL ||
             strcmp(info.status,
                    "indeterminate (inaccurate - reached max_iters)") != 0) {
    fail = "test_indeterminate_at_limit: max_iters 2: wrong status string";
  } else if (info.iter != 2) {
    fail = "test_indeterminate_at_limit: max_iters 2: expected 2 iterations";
  } else if (!inaccurate_sol_all_nan(sol, d->n, d->m)) {
    fail = "test_indeterminate_at_limit: max_iters 2: x, y, s not all NaN";
  } else if (!isnan(info.pobj) || !isnan(info.dobj) || !isnan(info.gap) ||
             !isnan(info.res_pri) || !isnan(info.res_dual) ||
             !isnan(info.res_infeas) || !isnan(info.res_unbdd_a) ||
             !isnan(info.res_unbdd_p) || !isnan(info.comp_slack)) {
    fail = "test_indeterminate_at_limit: max_iters 2: objective or residual "
           "fields not NaN";
  }

  if (!fail) {
    stgs->max_iters = 3;
    memset(&info, 0, sizeof(info));
    exitflag = scs(d, k, stgs, sol, &info);
    if (exitflag != SCS_SOLVED_INACCURATE) {
      fail = "test_indeterminate_at_limit: max_iters 3: expected "
             "SCS_SOLVED_INACCURATE";
    }
  }
  if (!fail) {
    stgs->max_iters = MAX_ITERS;
    exitflag = scs(d, k, stgs, sol, &info);
    if (exitflag != SCS_SOLVED) {
      fail = "test_indeterminate_at_limit: default max_iters: expected "
             "SCS_SOLVED";
    }
  }

  /* SCS_INDETERMINATE, then a warm-started solve on the same workspace */
  if (!fail) {
    stgs->max_iters = 2;
    w = scs_init(d, k, stgs);
    if (!w) {
      fail = "test_indeterminate_at_limit: scs_init failed";
    }
  }
  if (!fail) {
    exitflag = scs_solve(w, sol, &info, 0);
    if (exitflag != SCS_INDETERMINATE) {
      fail = "test_indeterminate_at_limit: workspace, max_iters 2: expected "
             "SCS_INDETERMINATE";
    }
  }
  if (!fail) {
    stgs->max_iters = MAX_ITERS;
    if (scs_update_settings(w, stgs) != 0) {
      fail = "test_indeterminate_at_limit: scs_update_settings failed";
    }
  }
  if (!fail) {
    exitflag = scs_solve(w, sol, &info, 1);
    if (exitflag != SCS_SOLVED || !(fabs(sol->x[0] + 2.0) < 1e-2)) {
      fail = "test_indeterminate_at_limit: warm start from the "
             "indeterminate result: expected SCS_SOLVED at x = -2";
    }
  }
  if (fail) {
    scs_printf("test_indeterminate_at_limit: status %li '%.*s', iter %li\n",
               (long)exitflag, (int)sizeof(info.status) - 1, info.status,
               (long)info.iter);
  }

  scs_finish(w);
  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

/*
 * SCS_FAILED for an iterate that is not finite at the limit: the box LP
 * warm-started at x = y = s = half the largest finite scs_float. The first
 * linear solve overflows, and the iterate is NaN from then on. This tests an
 * overflow, which every build reaches; a linear solver or BLAS that returns
 * a NaN takes the same path. Before this change max_iters = 1 returned
 * SCS_SOLVED_INACCURATE with a non-finite solution (tau is fixed at 1 in
 * iteration 0) and max_iters = 2 returned SCS_FAILED with the status
 * " (inaccurate - reached max_iters)". A linear solver that rejects the
 * overflowed system itself stops the solve with SCS_FAILED and "failure"
 * too, at another iteration count, so info.iter is not checked. The finite
 * runs of test_indeterminate_at_limit are the other side.
 */
static const char *test_nonfinite_at_limit(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag, max_iters, i;
#ifdef SFLOAT
  const scs_float big = FLT_MAX / 2;
#else
  const scs_float big = DBL_MAX / 2;
#endif
  const char *fail = SCS_NULL;

  inaccurate_box_lp(d, k);
  scs_set_default_settings(stgs);
  stgs->verbose = 1;
  stgs->warm_start = 1;
  sol->x = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  sol->y = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  sol->s = (scs_float *)scs_calloc(d->m, sizeof(scs_float));

  for (max_iters = 1; max_iters <= 2 && !fail; ++max_iters) {
    for (i = 0; i < d->n; ++i) {
      sol->x[i] = big;
    }
    for (i = 0; i < d->m; ++i) {
      sol->y[i] = sol->s[i] = big;
    }
    stgs->max_iters = max_iters;
    exitflag = scs(d, k, stgs, sol, &info);
    if (exitflag != SCS_FAILED || info.status_val != SCS_FAILED) {
      fail = "test_nonfinite_at_limit: expected SCS_FAILED";
    } else if (strcmp(info.status, "failure") != 0) {
      fail = "test_nonfinite_at_limit: expected the status \"failure\"";
    } else if (!inaccurate_sol_all_nan(sol, d->n, d->m)) {
      fail = "test_nonfinite_at_limit: x, y, s not all NaN";
    }
    if (fail) {
      scs_printf("test_nonfinite_at_limit: max_iters %li: status %li '%s'\n",
                 (long)max_iters, (long)exitflag, info.status);
    }
  }

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

#undef _INACCURATE_LP_VARS
#undef _INACCURATE_INFEAS_VARS
#undef _INACCURATE_UNBDD_VARS
