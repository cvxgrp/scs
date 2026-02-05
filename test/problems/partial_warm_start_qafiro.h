#include <string.h>

#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

/* Partial warm start tests on the qafiro QP (n=32, m=60).
 * Tests various partial initialization patterns.
 */
static const char *partial_warm_start_qafiro(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_int cold_iters;
  scs_int i;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  /* saved solution for partial warm starts */
  scs_float saved_x[32];
  scs_float saved_y[60];
  scs_float saved_s[60];

  /* data (same as qafiro_tiny_qp) */
  scs_float Ax[] = {
      -1.,    -1.06, -1.,    -0.301, -1.,    1.,     1.,     -1.,   1.,
      1.,     -1.,   1.,     -1.,    -1.,    -1.,    -1.06,  -1.,   -0.301,
      -1.,    -1.,   -1.06,  -1.,    -0.313, -1.,    -1.,    -0.96, -1.,
      -0.313, -1.,   -1.,    -0.86,  -1.,    -0.326, -1.,    1.,    -2.364,
      -1.,    1.,    -2.386, -1.,    1.,     -2.408, -1.,    1.,    -2.429,
      -1.,    1.,    -1.4,   -1.,    1.,     1.,     -1.,    1.,    -1.,
      -1.,    -0.43, -1.,    -1.,    -0.109, -1.,    1.,     1.,    -1.,
      1.,     1.,    -1.,    1.,     1.,     -1.,    1.,     -1.,   -1.,
      1.,     -0.43, -1.,    -0.109, -1.,    1.,     -0.43,  -1.,   -0.108,
      -1.,    1.,    -0.39,  -1.,    -0.108, -1.,    1.,     -0.37, -1.,
      -0.107, -1.,   1.,     -2.191, -1.,    1.,     -2.219, -1.,   1.,
      -2.249, -1.,   1.,     -2.279, -1.,    -1.,    -1.4,   -1.,   1.,
      1.,     -1.,   1.,     -1.,    -1.,    1.,     -1.};
  scs_int Ai[] = {0,  1,  16, 24, 28, 0,  15, 29, 0,  22, 30, 1,  26, 31, 4,
                  5,  12, 25, 32, 4,  5,  11, 25, 33, 4,  5,  9,  25, 34, 4,
                  5,  10, 25, 35, 12, 21, 36, 11, 21, 37, 9,  21, 38, 10, 21,
                  39, 4,  15, 40, 4,  23, 41, 5,  27, 42, 6,  7,  13, 22, 43,
                  7,  14, 44, 7,  24, 45, 7,  21, 46, 6,  26, 47, 2,  3,  17,
                  23, 48, 2,  3,  18, 23, 49, 2,  3,  19, 23, 50, 2,  3,  20,
                  23, 51, 17, 21, 52, 18, 21, 53, 19, 21, 54, 20, 21, 55, 2,
                  14, 56, 2,  25, 57, 3,  27, 58, 2,  59};
  scs_int Ap[] = {0,  5,  8,  11, 14, 19,  24,  29,  34,  37,  40,
                  43, 46, 49, 52, 55, 60,  63,  66,  69,  72,  77,
                  82, 87, 92, 95, 98, 101, 104, 107, 110, 113, 115};

  scs_float Px[] = {10., 1., 10., 1., 1., 10.};
  scs_int Pi[] = {0, 0, 1, 0, 1, 2};
  scs_int Pp[] = {0, 1, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};

  scs_float b[] = {
      0.00000000e+00, 0.00000000e+00, 4.40000000e+01, -2.22044605e-16,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00};
  scs_float c_data[] = {0., -0.4,  0., 0., 0., 0.,   0.,    0., 0., 0., 0.,
                        0., -0.32, 0., 0., 0., -0.6, 0.,    0., 0., 0., 0.,
                        0., 0.,    0., 0., 0., 0.,   -0.48, 0., 0., 10.};

  scs_int m_val = 60;
  scs_int n_val = 32;

  scs_float bl[] = {
      -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20,
      -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20,
      -1e+20, 0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      0.0,    0.0,    0.0,    0.0,    0.0,    0.0};
  scs_float bu[] = {0.0,   0.0,   0.0,   80.0,  500.0, 0.0,   0.0,   80.0,
                    500.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                    0.0,   310.0, 300.0, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20};
  scs_int bsize = 52;
  scs_int z = 8;

  scs_float opt = -1.5907818;
  /* end data */

  d->m = m_val;
  d->n = n_val;
  d->b = b;
  d->c = c_data;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));

  d->A->m = m_val;
  d->A->n = n_val;

  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  d->P->m = n_val;
  d->P->n = n_val;

  d->P->x = Px;
  d->P->i = Pi;
  d->P->p = Pp;

  k->bsize = bsize;
  k->bl = bl;
  k->bu = bu;
  k->z = z;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;
  stgs->eps_infeas = 0.;

  /* Step 1: Cold solve */
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: cold solve failed", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) {
    goto cleanup;
  }

  cold_iters = info.iter;
  scs_printf("partial_warm_start_qafiro: cold solve took %li iters\n",
             (long)cold_iters);

  /* Save solution */
  memcpy(saved_x, sol->x, n_val * sizeof(scs_float));
  memcpy(saved_y, sol->y, m_val * sizeof(scs_float));
  memcpy(saved_s, sol->s, m_val * sizeof(scs_float));

  /* Step 2: Keep x, NaN out s and y */
  memcpy(sol->x, saved_x, n_val * sizeof(scs_float));
  for (i = 0; i < m_val; ++i) {
    sol->s[i] = NAN;
    sol->y[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: x-only warm start failed", success);
  scs_printf(
      "partial_warm_start_qafiro: x-only warm start took %li iters\n",
      (long)info.iter);

  /* Step 3: Keep first half of x, NaN rest of x, s, y */
  for (i = 0; i < n_val; ++i) {
    sol->x[i] = (i < n_val / 2) ? saved_x[i] : NAN;
  }
  for (i = 0; i < m_val; ++i) {
    sol->s[i] = NAN;
    sol->y[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: half-x warm start failed", success);
  scs_printf(
      "partial_warm_start_qafiro: half-x warm start took %li iters\n",
      (long)info.iter);

  /* Step 4: Keep x and y, NaN out only s
   * s = b - Ax from optimal x gives optimal s, so this should help a lot
   */
  memcpy(sol->x, saved_x, n_val * sizeof(scs_float));
  memcpy(sol->y, saved_y, m_val * sizeof(scs_float));
  for (i = 0; i < m_val; ++i) {
    sol->s[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: x+y warm start (NaN s) failed",
            success);
  scs_printf(
      "partial_warm_start_qafiro: x+y warm start (NaN s) took %li iters\n",
      (long)info.iter);
  mu_assert(
      "partial_warm_start_qafiro: x+y (NaN s) should beat cold start",
      info.iter < cold_iters);

  /* Step 5: Keep x and s, NaN out only y */
  memcpy(sol->x, saved_x, n_val * sizeof(scs_float));
  memcpy(sol->s, saved_s, m_val * sizeof(scs_float));
  for (i = 0; i < m_val; ++i) {
    sol->y[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: x+s warm start (NaN y) failed",
            success);
  scs_printf(
      "partial_warm_start_qafiro: x+s warm start (NaN y) took %li iters\n",
      (long)info.iter);
  /* x+s with NaN y can be worse than cold start on some platforms,
   * so just check convergence here (no iteration count assertion). */

  /* Step 6: Keep all of x and y except first 2 entries each, NaN out s */
  memcpy(sol->x, saved_x, n_val * sizeof(scs_float));
  sol->x[0] = NAN;
  sol->x[1] = NAN;
  memcpy(sol->y, saved_y, m_val * sizeof(scs_float));
  sol->y[0] = NAN;
  sol->y[1] = NAN;
  for (i = 0; i < m_val; ++i) {
    sol->s[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: most x+y (NaN s) failed", success);
  scs_printf(
      "partial_warm_start_qafiro: most x+y (first 2 NaN, NaN s) took %li iters\n",
      (long)info.iter);
  if (strstr(info.lin_sys_solver, "indirect")) {
    mu_assert(
        "partial_warm_start_qafiro: most x+y (NaN s) should converge in <= cold iters",
        info.iter <= cold_iters);
  } else {
    mu_assert(
        "partial_warm_start_qafiro: most x+y (NaN s) should beat cold start",
        info.iter < cold_iters);
  }

  /* Step 7: Perturbed x and y, NaN s (simulates nearby problem solution) */
  for (i = 0; i < n_val; ++i) {
    sol->x[i] = saved_x[i] * 1.01; /* 1% perturbation */
  }
  for (i = 0; i < m_val; ++i) {
    sol->y[i] = saved_y[i] * 1.01;
    sol->s[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: perturbed x+y (1%%) failed", success);
  scs_printf(
      "partial_warm_start_qafiro: perturbed x+y 1%% (NaN s) took %li iters\n",
      (long)info.iter);
  if (strstr(info.lin_sys_solver, "indirect")) {
    mu_assert(
        "partial_warm_start_qafiro: perturbed x+y (1%%) should converge in <= cold iters",
        info.iter <= cold_iters);
  } else {
    mu_assert(
        "partial_warm_start_qafiro: perturbed x+y (1%%) should beat cold start",
        info.iter < cold_iters);
  }

  /* Step 8: Larger perturbation */
  for (i = 0; i < n_val; ++i) {
    sol->x[i] = saved_x[i] * 1.1; /* 10% perturbation */
  }
  for (i = 0; i < m_val; ++i) {
    sol->y[i] = saved_y[i] * 1.1;
    sol->s[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;
  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start_qafiro: perturbed x+y (10%%) failed", success);
  scs_printf(
      "partial_warm_start_qafiro: perturbed x+y 10%% (NaN s) took %li iters\n",
      (long)info.iter);
  if (strstr(info.lin_sys_solver, "indirect")) {
    mu_assert(
        "partial_warm_start_qafiro: perturbed x+y (10%%) should converge in <= cold iters",
        info.iter <= cold_iters);
  } else {
    mu_assert(
        "partial_warm_start_qafiro: perturbed x+y (10%%) should beat cold start",
        info.iter < cold_iters);
  }

cleanup:
  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
