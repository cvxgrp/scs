#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <string.h>

/*
 * In double precision the per-iteration refinement turns off at a check
 * whose plain solve is at rounding level, and does not turn on there
 * (REFINE_ROUNDING_EPS in glbopts.h), here on LPs where one of SCS's
 * residuals goes to zero (a certificate) and with it the tol that
 * refinement compares with.
 *
 * Each LP (n variables, 2 n rows, every row a nonnegative cone) is solved
 * for RC_ITERS iterations with eps_infeas 0, so that SCS keeps iterating on
 * the certificate (as test_inaccurate does), at a fixed scale
 * (adaptive_scale off, R_y = 1 / scale), and the test requires refinement
 * to be on for at most a tenth of the iterations. With the absolute floor
 * before 3.4.0 (tol at least 1e-12, whatever the size of the terms of the
 * residual) refinement turned on once the residual reached zero and stayed
 * on to the end, chasing rounding:
 *
 * - An unbounded LP (RC_N, RC_SEED, RC_SCALE; every linear solver): on for
 *   799 of the 1000 iterations with QDLDL (dense, Accelerate, cuDSS and CG
 *   did not turn on, before or after). Now 0.
 * - An infeasible LP (RC_CG_N, RC_CG_SEED, RC_CG_SCALE; the CG solvers): on
 *   for 703 iterations. Now 0. The direct solvers skip it: there K is badly
 *   scaled and their plain solves are well above rounding level, so
 *   refinement is rightly on.
 *
 * Skipped in single precision, where refinement has no rounding test and
 * stays on as before (on these LPs 949 and 800 iterations): there it also
 * lowers an error of the solution that the residual does not show, and
 * turning it off at rounding level lost verified solves (see
 * REFINE_ROUNDING_EPS). Measured with clang on Apple silicon (QDLDL, dense,
 * Accelerate, CPU CG) and with GCC on ARM (QDLDL, cuDSS, CPU CG).
 */

/* RC_N, RC_SEED, RC_SCALE and RC_CG_N, RC_CG_SEED, RC_CG_SCALE: the two
 * instances (see above). RC_ITERS: the iterations of each run: long enough
 * that refinement left on chasing rounding (70 to 80% of them before) is far
 * from the bar. */
#define RC_N (20)
#define RC_SEED (2)
#define RC_SCALE (1e-2)
#define RC_CG_N (30)
#define RC_CG_SEED (1)
#define RC_CG_SCALE (1e-4)
#define RC_ITERS (1000)

#ifndef SFLOAT
/* a random LP with A dense (m = 2 n rows), every row a nonnegative cone,
 * infeasible (y >= 0 with A'y = 0 and b'y < 0) or unbounded (a direction
 * dx with A dx <= 0 and c'dx < 0, and b with A x + s = b feasible: x = 0,
 * s = b > 0) */
static void rc_data(scs_int infeas, scs_int n, scs_int seed, ScsData *d,
                    ScsCone *k) {
  scs_int i, j, m = 2 * n;
  scs_float *v, vv = 0., t;
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
  k->l = m;
  ran_start(seed);
  for (j = 0; j < n; ++j) {
    A->p[j + 1] = (j + 1) * m;
    for (i = 0; i < m; ++i) {
      A->i[j * m + i] = i;
      A->x[j * m + i] = rand_scs_float();
    }
  }
  if (infeas) {
    /* y = v >= 0: remove A's columns' components along v, then b'v = -1 */
    v = (scs_float *)scs_calloc(m, sizeof(scs_float));
    for (i = 0; i < m; ++i) {
      v[i] = 0.5 + 0.5 * rand_scs_float();
      vv += v[i] * v[i];
    }
    for (j = 0; j < n; ++j) {
      t = 0.;
      for (i = 0; i < m; ++i) {
        t += v[i] * A->x[j * m + i];
      }
      for (i = 0; i < m; ++i) {
        A->x[j * m + i] -= t * v[i] / vv;
      }
    }
    t = 0.;
    for (i = 0; i < m; ++i) {
      d->b[i] = rand_scs_float();
      t += v[i] * d->b[i];
    }
    for (i = 0; i < m; ++i) {
      d->b[i] -= (t + 1.) * v[i] / vv;
    }
    for (j = 0; j < n; ++j) {
      d->c[j] = rand_scs_float();
    }
  } else {
    /* dx = v: each row of A gets A_i dx in [-1, -0.5], then c'dx = -1 */
    v = (scs_float *)scs_calloc(n, sizeof(scs_float));
    for (j = 0; j < n; ++j) {
      v[j] = rand_scs_float();
      vv += v[j] * v[j];
    }
    for (i = 0; i < m; ++i) {
      scs_float ad = 0., target = -(0.5 + 0.5 * rand_scs_float());
      for (j = 0; j < n; ++j) {
        ad += A->x[j * m + i] * v[j];
      }
      for (j = 0; j < n; ++j) {
        A->x[j * m + i] += (target - ad) * v[j] / vv;
      }
    }
    t = 0.;
    for (j = 0; j < n; ++j) {
      d->c[j] = rand_scs_float();
      t += d->c[j] * v[j];
    }
    for (j = 0; j < n; ++j) {
      d->c[j] -= (t + 1.) * v[j] / vv;
    }
    for (i = 0; i < m; ++i) {
      d->b[i] = 1. + rand_scs_float();
    }
  }
  scs_free(v);
}

/* solves one instance for RC_ITERS iterations at a fixed scale; returns an
 * error message or SCS_NULL */
static const char *rc_solve(const char *name, scs_int infeas, scs_int n,
                            scs_int seed, scs_float scale) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  ScsWork *w;
  scs_int status = SCS_FAILED, on_iters = -1, checks = -1, on = -1;
  rc_data(infeas, n, seed, d, k);
  scs_set_default_settings(stgs);
  stgs->scale = scale;
  stgs->adaptive_scale = 0;
  stgs->eps_infeas = 0.; /* keep iterating on the certificate */
  stgs->max_iters = RC_ITERS;
  stgs->verbose = 0;
  w = scs_init(d, k, stgs);
  if (w) {
    status = scs_solve(w, sol, &info, 0);
    on_iters = w->refine_on_iters;
    checks = w->refine_checks;
    on = w->refine_on;
    scs_finish(w);
  }
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  scs_printf("test_refine_certificate: %s, %s: status %li, %li iterations, "
             "%li checks, %li with refinement on, on at the end: %li\n",
             scs_get_lin_sys_method(), name, (long)status, (long)info.iter,
             (long)checks, (long)on_iters, (long)on);
  mu_assert("test_refine_certificate: did not run RC_ITERS iterations",
            info.iter == RC_ITERS);
  mu_assert("test_refine_certificate: refinement stayed on chasing rounding",
            on_iters >= 0 && 10 * on_iters <= RC_ITERS);
  return SCS_NULL;
}

#endif

static const char *test_refine_certificate(void) {
#ifdef SFLOAT
  return mu_skipped; /* no rounding test in single precision (see above) */
#else
  const char *fail = rc_solve("unbounded LP", 0, RC_N, RC_SEED, RC_SCALE);
  if (!fail && strstr(scs_get_lin_sys_method(), "indirect")) {
    fail = rc_solve("infeasible LP", 1, RC_CG_N, RC_CG_SEED, RC_CG_SCALE);
  }
  return fail;
#endif
}
