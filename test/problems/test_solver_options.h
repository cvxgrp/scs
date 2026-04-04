#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Shared setup for solver-options tests.
 * Problem: min x  s.t.  x >= 2  (n=1, m=1, l=1)
 * Optimal obj = 2.
 */
#define _OPTS_SETUP(bval, cval)                                                \
  do {                                                                         \
    scs_float _Ax[] = {-1.0};                                                  \
    scs_int _Ai[] = {0};                                                       \
    scs_int _Ap[] = {0, 1};                                                    \
    k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));                             \
    d = (ScsData *)scs_calloc(1, sizeof(ScsData));                             \
    stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));                  \
    sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));                   \
    d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));                      \
    d->b = (scs_float *)scs_calloc(1, sizeof(scs_float));                      \
    d->c = (scs_float *)scs_calloc(1, sizeof(scs_float));                      \
    d->m = 1; d->n = 1;                                                        \
    d->b[0] = (bval); d->c[0] = (cval);                                        \
    d->A->m = 1; d->A->n = 1;                                                  \
    d->A->x = (scs_float *)scs_calloc(1, sizeof(scs_float));                   \
    d->A->i = (scs_int *)scs_calloc(1, sizeof(scs_int));                       \
    d->A->p = (scs_int *)scs_calloc(2, sizeof(scs_int));                       \
    d->A->x[0] = _Ax[0]; d->A->i[0] = _Ai[0];                                 \
    d->A->p[0] = _Ap[0]; d->A->p[1] = _Ap[1];                                 \
    k->l = 1;                                                                  \
    scs_set_default_settings(stgs);                                            \
    stgs->eps_abs = 1e-6;                                                      \
    stgs->eps_rel = 1e-6;                                                      \
    stgs->verbose = 0;                                                         \
  } while (0)

#define _OPTS_CLEANUP()                                                        \
  do {                                                                         \
    SCS(free_sol)(sol);                                                        \
    scs_free(d->A->x);                                                         \
    scs_free(d->A->i);                                                         \
    scs_free(d->A->p);                                                         \
    scs_free(d->A);                                                            \
    scs_free(d->b);                                                            \
    scs_free(d->c);                                                            \
    scs_free(d);                                                               \
    scs_free(k);                                                               \
    scs_free(stgs);                                                            \
  } while (0)

/*
 * Test adaptive_scale=1: verify the solver converges to the correct answer.
 * Default settings already have adaptive_scale=1, but most problem tests
 * override it to 0 to get reproducible iteration counts.  This test uses a
 * random LP large enough (m=60, n=20) for the scale update heuristic to fire
 * at least once during the solve.
 */
static const char *test_adaptive_scale(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr;
  scs_int n = 20, m = 60;
  scs_int col_nnz = 4, nnz = n * col_nnz;
  const char *fail;

  k->z = 6;
  k->l = m - k->z;
  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, 4321);

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-5;
  stgs->eps_rel = 1e-5;
  stgs->adaptive_scale = 1;
  stgs->verbose = 0;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - SCS(dot)(d->c, opt_sol->x, d->n);

  mu_assert("test_adaptive_scale: expected SCS_SOLVED", exitflag == SCS_SOLVED);
  mu_assert("test_adaptive_scale: primal obj error too large", ABS(perr) < 1e-3);
  mu_assert("test_adaptive_scale: expected scale updates > 0",
            info.scale_updates > 0);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  return fail;
}

/*
 * Test acceleration_lookback=0: disables Anderson acceleration (w->accel=NULL
 * in init_work).  Verify the solver still converges to the correct answer.
 */
static const char *test_no_acceleration(void) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  _OPTS_SETUP(-2.0, 1.0);
  stgs->acceleration_lookback = 0;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_no_acceleration: expected SCS_SOLVED", exitflag == SCS_SOLVED);
  mu_assert("test_no_acceleration: no AA steps should be accepted",
            info.accepted_accel_steps == 0);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  _OPTS_CLEANUP();
  return fail;
}

/*
 * Test Type-II Anderson acceleration (acceleration_lookback < 0): the
 * negative lookback is the internal signal to use AA_REGULARIZATION_TYPE_2.
 * Verify convergence on a simple LP.
 */
static const char *test_type2_acceleration(void) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  _OPTS_SETUP(-2.0, 1.0);
  stgs->acceleration_lookback = -10;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_type2_acceleration: expected SCS_SOLVED",
            exitflag == SCS_SOLVED);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  _OPTS_CLEANUP();
  return fail;
}

/*
 * Test normalize=0 vs normalize=1: both should converge to the same primal
 * objective on the same LP (up to solver tolerance).
 */
static const char *test_normalize_off(void) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info0 = {0}, info1 = {0};
  scs_int exitflag;
  scs_float pobj0, pobj1;
  const char *fail;

  /* solve with normalize=1 */
  _OPTS_SETUP(-2.0, 1.0);
  stgs->normalize = 1;
  exitflag = scs(d, k, stgs, sol, &info1);
  mu_assert("test_normalize_off: normalize=1 failed", exitflag == SCS_SOLVED);
  pobj1 = info1.pobj;
  fail = verify_solution_correct(d, k, stgs, &info1, sol, exitflag);
  _OPTS_CLEANUP();
  if (fail) return fail;

  /* solve with normalize=0 */
  _OPTS_SETUP(-2.0, 1.0);
  stgs->normalize = 0;
  exitflag = scs(d, k, stgs, sol, &info0);
  mu_assert("test_normalize_off: normalize=0 failed", exitflag == SCS_SOLVED);
  pobj0 = info0.pobj;
  fail = verify_solution_correct(d, k, stgs, &info0, sol, exitflag);
  _OPTS_CLEANUP();
  if (fail) return fail;

  mu_assert("test_normalize_off: objectives differ between normalize=0 and 1",
            ABS(pobj0 - pobj1) < 1e-4);
  return 0;
}

#undef _OPTS_SETUP
#undef _OPTS_CLEANUP
