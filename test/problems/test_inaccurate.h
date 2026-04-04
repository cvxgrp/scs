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
  stgs->verbose = 0;

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
  stgs->max_iters = 500;
  stgs->verbose = 0;

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
  stgs->max_iters = 500;
  stgs->verbose = 0;

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

#undef _INACCURATE_LP_VARS
#undef _INACCURATE_INFEAS_VARS
#undef _INACCURATE_UNBDD_VARS
