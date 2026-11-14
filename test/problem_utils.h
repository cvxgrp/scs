#ifndef PUTILS_H_GUARD
#define PUTILS_H_GUARD

#include "cones.h"
#include "linalg.h"
#include "linsys.h"
#include "minunit.h"
#include "rng.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

#include <string.h>

#define _MAX_RAND_VAL (1073741823) /* 2^30 - 1 */

/* The tolerances and checks in this suite assume double precision: many sit
 * below what a single-precision solve can reach (machine epsilon is 1.2e-7, and
 * the KKT solve loses a further factor of 1 / rho_x). So in SFLOAT builds every
 * solve is run with tolerances no tighter than these, which keeps each test's
 * problem and expected status while asking for an achievable accuracy. A
 * tolerance of zero or less is left alone since it can be deliberate (disabling
 * a check) or the invalid input a validation test wants rejected. */
#ifdef SFLOAT
#define SFLOAT_TEST_EPS (1e-4)
#define SFLOAT_TEST_EPS_INFEAS (1e-4)
static void raise_tols_for_single_precision(ScsSettings *stgs) {
  if (!stgs) {
    return;
  }
  if (stgs->eps_abs > 0 && stgs->eps_abs < SFLOAT_TEST_EPS) {
    stgs->eps_abs = SFLOAT_TEST_EPS;
  }
  if (stgs->eps_rel > 0 && stgs->eps_rel < SFLOAT_TEST_EPS) {
    stgs->eps_rel = SFLOAT_TEST_EPS;
  }
  if (stgs->eps_infeas > 0 && stgs->eps_infeas < SFLOAT_TEST_EPS_INFEAS) {
    stgs->eps_infeas = SFLOAT_TEST_EPS_INFEAS;
  }
}
#define scs_init(d, k, stgs)                                                   \
  (raise_tols_for_single_precision(stgs), scs_init(d, k, stgs))
#define scs(d, k, stgs, sol, info)                                             \
  (raise_tols_for_single_precision(stgs), scs(d, k, stgs, sol, info))
/* double, single: recomputing a residual or objective outside the solver
 * disagrees by roundoff, about 1e-7 of the magnitudes involved */
#define VERIFY_TOL(dbl, sgl) (sgl)
#else
#define VERIFY_TOL(dbl, sgl) (dbl)
#endif
/* how close a solve must get to a known optimal objective or solution */
#define SOL_TOL VERIFY_TOL(1e-4, 1e-2)
/* the solver stops as soon as it meets its tolerances, so its roundoff can put
 * a recomputed residual just outside them */
#define FEAS_SLACK VERIFY_TOL(1., 1.1)

/* uniform random number in [-1,1] */
static scs_float rand_scs_float(void) {
  return 2 * (((scs_float)ran_arr_next()) / _MAX_RAND_VAL) - 1; /* in [-1, 1] */
}

void gen_random_prob_data(scs_int nnz, scs_int col_nnz, ScsData *d, ScsCone *k,
                          ScsSolution *opt_sol, scs_int seed) {
  scs_int n = d->n;
  scs_int m = d->m;
  ScsMatrix *A = d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  scs_float *b = d->b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *c = d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *x = opt_sol->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *y = opt_sol->y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *s = opt_sol->s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  /* temporary variables */
  scs_float *z = (scs_float *)scs_calloc(m, sizeof(scs_float));
  ScsConeWork *tmp_cone_work;
  scs_int i, j, r, rn, rm;

  A->i = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  A->p = (scs_int *)scs_calloc((n + 1), sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(nnz, sizeof(scs_float));
  A->n = d->n;
  A->m = d->m;
  /* seed before the first draw, so that the generated instance is a pure
   * function of (seed, n, m, nnz, col_nnz, k) rather than of whatever state
   * the generator was left in by an earlier caller */
  ran_start(seed);

  /* y, s >= 0 and y'*s = 0 */
  for (i = 0; i < m; i++) {
    y[i] = z[i] = rand_scs_float();
  }
  tmp_cone_work = SCS(init_cone)(k, m);
  SCS(proj_dual_cone)(y, tmp_cone_work, SCS_NULL, SCS_NULL);
  SCS(finish_cone(tmp_cone_work));

  for (i = 0; i < m; i++) {
    b[i] = s[i] = y[i] - z[i];
  }

  for (i = 0; i < n; i++) {
    x[i] = rand_scs_float();
  }

  /*
   c = -A'*y
   b = A*x + s
   */
  A->p[0] = 0;
  for (j = 0; j < n; j++) { /* column */
    r = 0;
    for (i = 0; i < m && r < col_nnz; ++i) {
      /* generate a unique sorted array via Knuths alg */
      rn = m - i;
      rm = col_nnz - r;
      if ((ran_arr_next() % rn) < rm) {
        A->x[r + j * col_nnz] = rand_scs_float();
        A->i[r + j * col_nnz] = i;
        b[i] += A->x[r + j * col_nnz] * x[j];
        c[j] -= A->x[r + j * col_nnz] * y[i];
        r++;
      }
    }
    A->p[j + 1] = (j + 1) * col_nnz;
  }
  scs_free(z);
}

/* Tolerance on the cone distances that the tests assert.
 *
 * A distance is measured with the cone projection itself, so it is only as
 * accurate as that projection. Every projection is closed-form or converges
 * to rounding, except the logdet cone's: that one is an iterative solve that
 * falls back from Newton to an interior-point method, and the fallback can
 * stop with a primal residual near 1e-4 (log_cone_wrapper accepts residuals
 * up to 1e-2). A point on the boundary of the cone, inside it to 1e-16, has
 * measured 2e-5 away that way.
 *
 * Which solves end on such a point is decided by rounding. With the data of
 * test_spectral_metric_updates perturbed by 1e-13, about one solve in a
 * thousand measured 2e-5 and none measured more than 2.2e-5, so a limit of
 * 1e-5 fails whenever a change to the solver, the BLAS or the CPU moves the
 * iterates onto one. 1e-3 is about 50 times the largest value seen. */
static scs_float cone_dist_tol(const ScsCone *k) {
#ifdef USE_SPECTRAL_CONES
  if (k->dsize > 0) {
    return 1e-3;
  }
#endif
  (void)k;
  return 1e-5;
}

/* distance of v (length len) from the cone: in single precision a projection
 * is accurate to about 1e-7 of the vector's size, so the limit grows with it */
#define CONE_DIST_TOL(k, v, len)                                               \
  VERIFY_TOL(cone_dist_tol(k),                                                 \
             MAX(cone_dist_tol(k), 1e-4 * (1 + NORM(v, len))))

/* single precision: the rounding in recomputing M v, where v is a certificate
 * (scaled so that b'y or c'x is -1, so its size is set by the problem) and the
 * product nearly cancels: a few machine epsilons of ||M||_F ||v|| */
#define PRODUCT_ROUND(M, v, len)                                               \
  (1e-5 * (1 + NORM((M)->x, (M)->p[(M)->n]) * NORM(v, len)))
/* the same for |Ax + s| and |Px| of an unboundedness certificate x */
#define UNBDD_A_ROUND(d, x, s, n, m)                                           \
  (PRODUCT_ROUND((d)->A, x, n) + 1e-5 * NORM(s, m))
#define UNBDD_P_ROUND(d, x, n) ((d)->P ? PRODUCT_ROUND((d)->P, x, n) : 0.)

static scs_float get_dual_cone_dist(const scs_float *y, ScsConeWork *c,
                                    scs_int m) {
  scs_float dist;
  scs_float *t = (scs_float *)scs_calloc(m, sizeof(scs_float));
  memcpy(t, y, m * sizeof(scs_float));
  SCS(proj_dual_cone)(t, c, SCS_NULL, SCS_NULL);
  dist = SCS(norm_inf_diff)(t, y, m);
  scs_free(t);
  return dist;
}

/* via moreau */
static scs_float get_pri_cone_dist(const scs_float *s, ScsConeWork *c,
                                   scs_int m) {
  scs_float dist;
  scs_float *t = (scs_float *)scs_calloc(m, sizeof(scs_float));
  memcpy(t, s, m * sizeof(scs_float));
  SCS(scale_array)(t, -1.0, m);
  SCS(proj_dual_cone)(t, c, SCS_NULL, SCS_NULL);
  dist = SCS(norm_inf)(t, m); /* ||s - Pi_c(s)|| = ||Pi_c*(-s)|| */
  scs_free(t);
  return dist;
}

const char *verify_solution_correct(ScsData *d, ScsCone *k, ScsSettings *stgs,
                                    ScsInfo *info, ScsSolution *sol,
                                    scs_int status) {
  scs_int n = d->n, m = d->m;
  scs_float *x = sol->x;
  scs_float *y = sol->y;
  scs_float *s = sol->s;

  scs_float *c = d->c;
  scs_float *b = d->b;

  scs_float *primal = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *ax = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *dual = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *px = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *aty = (scs_float *)scs_calloc(n, sizeof(scs_float));

  scs_float res_pri, res_dual, res_infeas, res_unbdd_a, res_unbdd_p;
  scs_float ctx, bty, xt_p_x, gap, pobj, dobj, sty;
  scs_float grl, prl, drl;

  scs_float sdist = NAN, ydist = NAN;

  ScsConeWork *cone_work = SCS(init_cone)(k, m);

  /**************** PRIMAL *********************/
  memset(ax, 0, m * sizeof(scs_float));
  /* Ax */
  SCS(accum_by_a)(d->A, x, ax);

  memcpy(primal, ax, m * sizeof(scs_float));
  /* Ax + s */
  SCS(add_scaled_array)(primal, s, m, 1.);

  /* unbounded residual |Ax + s| */
  res_unbdd_a = NORM(primal, m);

  /* Ax + s - b */
  SCS(add_scaled_array)(primal, b, m, -1.0);

  res_pri = NORM(primal, m);

  /**************** DUAL *********************/
  memset(px, 0, n * sizeof(scs_float));
  if (d->P) {
    /* px = Px */
    SCS(accum_by_p)(d->P, x, px);
    xt_p_x = SCS(dot)(px, x, n);
    res_unbdd_p = NORM(px, n);
  } else {
    xt_p_x = 0;
    res_unbdd_p = 0;
  }

  memset(aty, 0, n * sizeof(scs_float));
  /* aty = A'y */
  SCS(accum_by_atrans)(d->A, y, aty);
  res_infeas = NORM(aty, n);

  memcpy(dual, aty, n * sizeof(scs_float));
  /* Px + A'y */
  SCS(add_scaled_array)(dual, px, n, 1.);
  /* Px + A'y + c */
  SCS(add_scaled_array)(dual, c, n, 1.0);

  res_dual = NORM(dual, n);

  /**************** CONES *****************/

  if (status == SCS_SOLVED || status == SCS_UNBOUNDED) {
    sdist = get_pri_cone_dist(sol->s, cone_work, m);
  }
  if (status == SCS_SOLVED || status == SCS_INFEASIBLE) {
    ydist = get_dual_cone_dist(sol->y, cone_work, m);
  }

  /**************** OTHERS *****************/
  sty = SCS(dot)(y, s, m);

  bty = SCS(dot)(y, b, m);
  ctx = SCS(dot)(x, c, n);

  gap = ABS(xt_p_x + ctx + bty);
  pobj = xt_p_x / 2. + ctx;
  dobj = -xt_p_x / 2. - bty;

  /************** OPTIMALITY ****************/

  /* TODO: the MAX expansion computes these norms many times */
  grl = MAX(MAX(ABS(xt_p_x), ABS(ctx)), ABS(bty));
  prl = MAX(MAX(NORM(b, m), NORM(s, m)), NORM(ax, m));
  drl = MAX(MAX(NORM(c, n), NORM(px, n)), NORM(aty, n));

  /**************** CLEANUP *****************/
  scs_free(primal);
  scs_free(dual);
  scs_free(px);
  scs_free(ax);
  scs_free(aty);
  SCS(finish_cone)(cone_work);

  /**************** ASSERTS *****************/
  if (status == SCS_SOLVED) {
    /* reported-vs-recomputed consistency: the recomputation here uses a
     * different BLAS summation order than the solver's internal pipeline,
     * so the discrepancy scales with the magnitude of the accumulated
     * quantities (~ the objective scale), like the gap and objective
     * checks below already do */
    mu_assert_less("Primal residual ERROR", ABS(res_pri - info->res_pri),
                   VERIFY_TOL(1e-10, 1e-4) * (1 + ABS(pobj)));
    mu_assert_less("Dual residual ERROR", ABS(res_dual - info->res_dual),
                   VERIFY_TOL(1e-10, 1e-4) * (1 + ABS(pobj)));
    mu_assert_less("Gap ERROR", ABS(gap - info->gap),
                   VERIFY_TOL(1e-7 * (1 + ABS(gap)), 1e-4 * (1 + grl)));
    mu_assert_less("Primal obj ERROR", ABS(pobj - info->pobj),
                   VERIFY_TOL(1e-9, 1e-4) * (1 + ABS(pobj)));
    mu_assert_less("Dual obj ERROR", ABS(dobj - info->dobj),
                   VERIFY_TOL(1e-9, 1e-4) * (1 + ABS(dobj)));
    /* slightly looser tol */
    mu_assert_less("Complementary slackness ERROR", ABS(sty),
                   VERIFY_TOL(5e-8, 1e-4) * MAX(NORM(s, m), NORM(y, m)));
    mu_assert_less("s cone dist ERROR", ABS(sdist), CONE_DIST_TOL(k, s, m));
    mu_assert_less("y cone dist ERROR", ABS(ydist), CONE_DIST_TOL(k, y, m));

    mu_assert_less("Primal feas ERROR", res_pri,
                   FEAS_SLACK * (stgs->eps_abs + stgs->eps_rel * prl));
    mu_assert_less("Dual feas ERROR", res_dual,
                   FEAS_SLACK * (stgs->eps_abs + stgs->eps_rel * drl));
    mu_assert_less("Gap feas ERROR", gap,
                   FEAS_SLACK * (stgs->eps_abs + stgs->eps_rel * grl));

  } else if (status == SCS_INFEASIBLE) {
    mu_assert_less("Infeas ERROR", ABS(res_infeas - info->res_infeas),
                   VERIFY_TOL(1e-8, PRODUCT_ROUND(d->A, y, m)));
    /* the solver scales y so that b'y = -1; recomputing the sum rounds in
     * proportion to its terms, which can be far larger than 1 */
    mu_assert_less("bty ERROR", ABS(bty + 1),
                   VERIFY_TOL(1e-12, 1e-5 * (1 + NORM(b, m) * NORM(y, m))));
    mu_assert_less("y cone dist ERROR", ABS(ydist), CONE_DIST_TOL(k, y, m));
    mu_assert_less("Infeas invalid ERROR", res_infeas,
                   FEAS_SLACK * stgs->eps_infeas +
                       VERIFY_TOL(0., PRODUCT_ROUND(d->A, y, m)));

  } else if (status == SCS_UNBOUNDED) {
    mu_assert_less("Unbdd_a ERROR", ABS(res_unbdd_a - info->res_unbdd_a),
                   VERIFY_TOL(1e-8, UNBDD_A_ROUND(d, x, s, n, m)));
    mu_assert_less("Unbdd_p ERROR", ABS(res_unbdd_p - info->res_unbdd_p),
                   VERIFY_TOL(1e-8, UNBDD_P_ROUND(d, x, n)));
    /* as for b'y above */
    mu_assert_less("ctx ERROR", ABS(ctx + 1),
                   VERIFY_TOL(1e-12, 1e-5 * (1 + NORM(c, n) * NORM(x, n))));
    mu_assert_less("s cone dist ERROR", ABS(sdist), CONE_DIST_TOL(k, s, m));
    mu_assert_less("Unbounded P invalid ERROR", res_unbdd_p,
                   FEAS_SLACK * stgs->eps_infeas +
                       VERIFY_TOL(0., UNBDD_P_ROUND(d, x, n)));
    mu_assert_less("Unbounded A invalid ERROR", res_unbdd_a,
                   FEAS_SLACK * stgs->eps_infeas +
                       VERIFY_TOL(0., UNBDD_A_ROUND(d, x, s, n, m)));

  } else {
    return "INVALID STATUS";
  }
  return 0;
}

#endif
