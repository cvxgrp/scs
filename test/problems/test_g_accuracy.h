#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * SCS solves K g = [c; -b], K = [R_x + P, A'; A, -R_y], once per
 * factorization and uses g in every iteration until the next one, so it
 * checks the residual of that solve and refines g until it is at most
 * G_SOLVE_TOL relative to the terms (glbopts.h). This test checks that
 * guarantee, for every linear solver the test binary is built with, on a
 * badly scaled QP that is not normalized (so the KKT matrix keeps the bad
 * scaling). The first check, on the g of the first factorization (a
 * one-iteration solve), is the one a plain solve fails (its residual is about
 * 1e-7). The second, after a solve with scale updates, guards the update
 * path: g recomputed with the new R and still accurate. A plain solve
 * usually passes it too, because the scale updates improve the conditioning.
 */

/* ||[c; -b] - K g||_inf over the largest inf-norm among [c; -b], R g, A' g_y,
 * P g_x + A' g_y and A g_x: the measure the solver checks against
 * G_SOLVE_TOL (g_residual in src/scs.c), written out again here */
static scs_float g_accuracy_rel_res(const ScsWork *w, const scs_float *g) {
  scs_int i, n = w->d->n, m = w->d->m;
  scs_float *ax = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *aty = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *px = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float res = 0., scale = 0., ri;
  SCS(accum_by_a)(w->d->A, g, ax);
  SCS(accum_by_atrans)(w->d->A, &(g[n]), aty);
  if (w->d->P) {
    SCS(accum_by_p)(w->d->P, g, px);
  }
  for (i = 0; i < n; ++i) {
    ri = w->d->c[i] - (w->diag_r[i] * g[i] + px[i] + aty[i]);
    res = MAX(res, ABS(ri));
    scale = MAX(scale, MAX(ABS(w->d->c[i]), ABS(w->diag_r[i] * g[i])));
    scale = MAX(scale, MAX(ABS(px[i] + aty[i]), ABS(aty[i])));
  }
  for (i = 0; i < m; ++i) {
    ri = -w->d->b[i] - (ax[i] - w->diag_r[n + i] * g[n + i]);
    res = MAX(res, ABS(ri));
    scale = MAX(scale, MAX(ABS(w->d->b[i]), ABS(ax[i])));
    scale = MAX(scale, ABS(w->diag_r[n + i] * g[n + i]));
  }
  scs_free(ax);
  scs_free(aty);
  scs_free(px);
  return res / MAX(scale, 1e-300);
}

/* residual of a plain solve of the same system, to show what the
 * refinement had to fix (printed only: it depends on the linear solver) */
static scs_float g_accuracy_unrefined(ScsWork *w) {
  scs_int i, n = w->d->n, m = w->d->m;
  scs_float *g0 = (scs_float *)scs_calloc(n + m, sizeof(scs_float)), res;
  memcpy(g0, w->d->c, n * sizeof(scs_float));
  for (i = 0; i < m; ++i) {
    g0[n + i] = -w->d->b[i];
  }
  scs_solve_lin_sys(w->p, g0, SCS_NULL, CG_BEST_TOL);
  res = g_accuracy_rel_res(w, g0);
  scs_free(g0);
  return res;
}

/* how badly the test data is scaled: single precision cannot solve the
 * double-precision instance to any accuracy, refined or not */
#ifdef SFLOAT
#define G_ACC_ROW_SCALE (1e2)
#define G_ACC_COL_SCALE (1e1)
#else
#define G_ACC_ROW_SCALE (1e4)
#define G_ACC_COL_SCALE (1e3)
#endif

/* the solver refines g until this measure is at most G_SOLVE_TOL; recomputing
 * it rounds differently, hence the factor of 10. A plain solve of this
 * problem measures about 1e-7 (QDLDL, Accelerate, CG, cuDSS) in double
 * precision, so the limit still tells a refined g from an unrefined one. */
#define G_ACC_LIMIT (10 * G_SOLVE_TOL)

static const char *test_g_accuracy(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  ScsWork *w;
  scs_int n = 80, m = 240, col_nnz = 9, nnz = n * col_nnz, i, j;
  scs_float res, res0;
  const char *fail = SCS_NULL;

  k->z = 24;
  k->l = m - k->z;
  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, 4321);
  /* badly scaled: every fifth row of A (and of b) times G_ACC_ROW_SCALE,
   * every seventh column (and entry of c) times G_ACC_COL_SCALE */
  for (j = 0; j < n; ++j) {
    for (i = d->A->p[j]; i < d->A->p[j + 1]; ++i) {
      if (d->A->i[i] % 5 == 0) {
        d->A->x[i] *= G_ACC_ROW_SCALE;
      }
      if (j % 7 == 0) {
        d->A->x[i] *= G_ACC_COL_SCALE;
      }
    }
    if (j % 7 == 0) {
      d->c[j] *= G_ACC_COL_SCALE;
    }
  }
  for (i = 0; i < m; ++i) {
    if (i % 5 == 0) {
      d->b[i] *= G_ACC_ROW_SCALE;
    }
  }
  /* diagonal P, so that the P term is checked too */
  d->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->P->m = d->P->n = n;
  d->P->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  d->P->i = (scs_int *)scs_calloc(n, sizeof(scs_int));
  d->P->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    d->P->p[j + 1] = j + 1;
    d->P->i[j] = j;
    d->P->x[j] = (j % 3 == 0) ? 1e2 : 1e-3;
  }

  scs_set_default_settings(stgs);
  stgs->normalize = 0;
  stgs->max_iters = 1;
  stgs->verbose = 0;

  w = scs_init(d, k, stgs);
  mu_assert("test_g_accuracy: scs_init returned NULL", w != SCS_NULL);
  /* a solve computes g at its start (update_work); one iteration makes no
   * scale update, so g is the one of the first factorization */
  scs_solve(w, sol, &info, 0);
  res = g_accuracy_rel_res(w, w->g);
  res0 = g_accuracy_unrefined(w);
  scs_printf("test_g_accuracy: first factorization, g relative residual "
             "%.2e (plain solve %.2e)\n",
             res, res0);
  if (info.scale_updates != 0) {
    fail = "test_g_accuracy: unexpected scale update in one iteration";
  } else if (res > G_ACC_LIMIT) {
    fail = "test_g_accuracy: g not accurate for the first factorization";
  }

  if (!fail) {
    /* scale updates refactorize and solve for g again with the new R */
    w->stgs->max_iters = 2000;
    scs_solve(w, sol, &info, 0);
    res = g_accuracy_rel_res(w, w->g);
    res0 = g_accuracy_unrefined(w);
    scs_printf("test_g_accuracy: after %li iterations and %li scale updates, "
               "g relative residual %.2e (plain solve %.2e)\n",
               (long)info.iter, (long)info.scale_updates, res, res0);
    if (info.scale_updates < 1) {
      fail = "test_g_accuracy: no scale update, the second check is vacuous";
    } else if (res > G_ACC_LIMIT) {
      fail = "test_g_accuracy: g not accurate after scale updates";
    }
  }

  scs_finish(w);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  return fail;
}
