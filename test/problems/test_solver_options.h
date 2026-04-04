#include <string.h>

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

/*
 * Test normalize roundtrip: verify that un_normalize_sol correctly inverts
 * the Ruiz equilibration by solving the same problem with normalize=0 and
 * normalize=1 and comparing the full primal solution vectors.
 *
 * The problem uses A entries with different magnitudes (1, 10, 100) so that
 * Ruiz scaling actually changes D/E significantly.
 *
 * Problem: min x1 + x2 + x3
 *          subject to  x1 >= 1
 *                     x2 >= 1
 *                     x3 >= 1
 * with A = diag(-1, -10, -100), b = [-1, -10, -100], c = [1, 1, 1].
 * Solution: x1 = x2 = x3 = 1, obj = 3.
 */
static const char *test_normalize_roundtrip(void) {
  ScsSettings *stgs0, *stgs1;
  ScsData *d0, *d1;
  ScsCone *k0, *k1;
  ScsSolution *sol0, *sol1;
  ScsInfo info0 = {0}, info1 = {0};
  scs_int exitflag0, exitflag1, i;
  const char *fail;

  scs_float Ax[] = {-1.0, -10.0, -100.0};
  scs_int Ai[]   = {0, 1, 2};
  scs_int Ap[]   = {0, 1, 2, 3};
  scs_float b[]  = {-1.0, -10.0, -100.0};
  scs_float c[]  = {1.0, 1.0, 1.0};
  scs_int m = 3, n = 3;

  /* helper macro to init data/settings/cone/sol */
#define _SETUP_NRM(d, k, stgs, sol) \
  do { \
    (d) = (ScsData *)scs_calloc(1, sizeof(ScsData)); \
    (d)->m = m; (d)->n = n; (d)->b = b; (d)->c = c; \
    (d)->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix)); \
    (d)->A->m = m; (d)->A->n = n; \
    (d)->A->x = Ax; (d)->A->i = Ai; (d)->A->p = Ap; \
    (k) = (ScsCone *)scs_calloc(1, sizeof(ScsCone)); \
    (k)->l = 3; \
    (stgs) = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings)); \
    scs_set_default_settings(stgs); \
    (stgs)->eps_abs = 1e-7; (stgs)->eps_rel = 1e-7; \
    (stgs)->verbose = 0; \
    (sol) = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution)); \
  } while (0)

#define _CLEANUP_NRM(d, k, stgs, sol) \
  do { \
    SCS(free_sol)(sol); \
    scs_free((d)->A); \
    scs_free(d); \
    scs_free(k); \
    scs_free(stgs); \
  } while (0)

  _SETUP_NRM(d1, k1, stgs1, sol1);
  stgs1->normalize = 1;
  exitflag1 = scs(d1, k1, stgs1, sol1, &info1);
  mu_assert("test_normalize_roundtrip: normalize=1 failed", exitflag1 == SCS_SOLVED);
  fail = verify_solution_correct(d1, k1, stgs1, &info1, sol1, exitflag1);
  if (fail) { _CLEANUP_NRM(d1, k1, stgs1, sol1); return fail; }

  _SETUP_NRM(d0, k0, stgs0, sol0);
  stgs0->normalize = 0;
  exitflag0 = scs(d0, k0, stgs0, sol0, &info0);
  mu_assert("test_normalize_roundtrip: normalize=0 failed", exitflag0 == SCS_SOLVED);
  fail = verify_solution_correct(d0, k0, stgs0, &info0, sol0, exitflag0);
  if (fail) { _CLEANUP_NRM(d0, k0, stgs0, sol0); return fail; }

  /* Compare primal solution vectors */
  for (i = 0; i < n; ++i) {
    mu_assert("test_normalize_roundtrip: x differs between normalize=0 and 1",
              ABS(sol0->x[i] - sol1->x[i]) < 1e-3);
  }
  mu_assert("test_normalize_roundtrip: objectives differ",
            ABS(info0.pobj - info1.pobj) < 1e-4);

  _CLEANUP_NRM(d0, k0, stgs0, sol0);
  _CLEANUP_NRM(d1, k1, stgs1, sol1);

#undef _SETUP_NRM
#undef _CLEANUP_NRM
  return 0;
}

/*
 * Test scs_version(): the public API function should return a non-NULL,
 * non-empty version string.
 */
static const char *test_scs_version(void) {
  const char *ver = scs_version();
  mu_assert("scs_version: should not return NULL", ver != SCS_NULL);
  mu_assert("scs_version: version string should be non-empty", ver[0] != '\0');
  return 0;
}

/*
 * Test time_limit_secs: setting an extremely small time limit (1e-10 s)
 * guarantees the timer fires before any iteration completes.  The solver must
 * exit with an inaccurate status code and append "(inaccurate - reached
 * time_limit_secs)" to info.status.
 */
static const char *test_time_limit_secs(void) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag;

  _OPTS_SETUP(-2.0, 1.0);
  stgs->time_limit_secs = 1e-10; /* fires before first iteration */
  stgs->max_iters = 100000;       /* large enough so max_iters is not the cause */

  exitflag = scs(d, k, stgs, sol, &info);

  /* Must return an inaccurate code (not SCS_FAILED) */
  mu_assert("test_time_limit_secs: should not return SCS_FAILED",
            exitflag != SCS_FAILED);
  mu_assert("test_time_limit_secs: expected inaccurate exit code",
            exitflag == SCS_SOLVED_INACCURATE ||
            exitflag == SCS_INFEASIBLE_INACCURATE ||
            exitflag == SCS_UNBOUNDED_INACCURATE);
  mu_assert("test_time_limit_secs: status should mention time_limit_secs",
            strstr(info.status, "time_limit_secs") != SCS_NULL);

  _OPTS_CLEANUP();
  return 0;
}

/*
 * Test explicit warm start: manually set sol->x, sol->y, sol->s to values
 * far from the optimum, then solve with warm_start=1.  The solver should
 * still converge to the correct answer.
 *
 * Uses the 3-step API: scs_init / scs_solve / scs_finish.
 *
 * Problem: min x  s.t.  x >= 2  (same LP as _OPTS_SETUP)
 * Optimal: x=2, obj=2.
 */
static const char *test_warm_start(void) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  _OPTS_SETUP(-2.0, 1.0);

  ScsWork *w = scs_init(d, k, stgs);
  mu_assert("test_warm_start: scs_init failed", w != SCS_NULL);

  /* Cold solve first to verify baseline */
  exitflag = scs_solve(w, sol, &info, 0);
  mu_assert("test_warm_start: cold solve failed", exitflag == SCS_SOLVED);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) { scs_finish(w); _OPTS_CLEANUP(); return fail; }

  /* Now set sol to a bad initial guess: x=100, y=0, s=0 */
  sol->x[0] = 100.0;
  sol->y[0] = 0.0;
  sol->s[0] = 0.0;

  /* Warm-start solve with the bad initial guess */
  exitflag = scs_solve(w, sol, &info, 1);
  mu_assert("test_warm_start: warm solve failed", exitflag == SCS_SOLVED);
  mu_assert("test_warm_start: wrong objective after warm start",
            ABS(info.pobj - 2.0) < 1e-4);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  scs_finish(w);
  _OPTS_CLEANUP();
  return fail;
}

#undef _OPTS_SETUP
#undef _OPTS_CLEANUP
