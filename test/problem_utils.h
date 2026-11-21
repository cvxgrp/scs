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

/* Tolerance on |s'y| of a solved problem.
 *
 * s and y come from one projection (Moreau), so s'y is zero in exact
 * arithmetic, and its rounding is that projection's: each part is off by some
 * unit roundoffs u of the vector projected, which is about as large as the
 * larger of ||s||_2 and ||y||_2. So |s'y| is a number of u times the square
 * of that, whatever the problem and its scaling: at most 40 u in single
 * precision and 60 u in double over the whole suite (Mac: QDLDL, CG, dense
 * and Accelerate, with and without FMA contraction; x86-64 with OpenBLAS:
 * -march=native and -march=x86-64 builds on Zen 3, Zen 4, Cascade Lake and
 * Granite Rapids).
 *
 * Single precision checks |s'y| against 1e-4 (the tolerance the harness
 * raises single-precision solves to) times max(||s||_2, ||y||_2)^2, about
 * 40 times the largest rounding seen. It was 1e-4 max(||s||_inf,
 * ||y||_inf), which grows like one of the vectors where s'y grows like both:
 * on test_complex_psd_multiplet_solve (||s||_2 = 10, ||y||_2 = 0.16) it was
 * 370 u of ||s||_2 ||y||_2, below the projection's rounding there (up to
 * 440 u on x86 hosts), and correct solves failed it.
 *
 * What it catches: an error of relative size e in s or y moves s'y by up to
 * e ||s||_2 ||y||_2 (projecting onto the second-order cone with t 0.1% too
 * large gave 1.6 times the tolerance on test_mixed_cones). Where it does not
 * work: when one of s and y is k times smaller than the other, an error in
 * the smaller one shows only above 1e-4 k of its size (complex_PSD, k = 27,
 * with the PSD projection's eigenvalues 1% too large: 3% of the tolerance;
 * the cone distance check caught it). Not ||s||_2 ||y||_2 instead: where s or
 * y is zero at the solution, its computed value is rounding alone and can lie
 * along the other (lp_update: ||s||_2 = 2e-6 and s'y = ||s||_2 ||y||_2). Not
 * inf-norms: |s'y| <= ||s||_2 ||y||_2, but only m ||s||_inf ||y||_inf.
 *
 * Double precision keeps 5e-8 max(||s||_inf, ||y||_inf), at least 7,000
 * times its rounding over the suite. */
static scs_float comp_slack_tol(const scs_float *s, const scs_float *y,
                                scs_int m) {
#ifdef SFLOAT
  scs_float sz = MAX(SCS(norm_2)(s, m), SCS(norm_2)(y, m));
  return 1e-4 * sz * sz;
#else
  return 5e-8 * MAX(NORM(s, m), NORM(y, m));
#endif
}

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
    mu_assert_less("Complementary slackness ERROR", ABS(sty),
                   comp_slack_tol(s, y, m));
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

/* Warm-start check: a warm start makes SCS start from the (x, y, s) it is
 * given, so the first iterate of a solve warm-started from a solution is close
 * to that solution, and the first iterate of a cold start is not. For each of
 * the blocks x, y and s this compares the two first iterates' distances from
 * sol, ||x1 - x||_inf and so on, and returns the largest ratio, warm over
 * cold. Both solves use stgs with max_iters = 1, so the result does not
 * depend on where a solve happens to stop; sol is not changed.
 *
 * A test asserts that the ratio is at most WARM_START_FIRST_ITER_TOL, set
 * between what a working warm start and a broken one give. Measured on
 * hs21_tiny_qp (Mac clang and Spark gcc, double and single precision, QDLDL,
 * CG, dense LAPACK and Accelerate, 2,000 starting scales around the default,
 * also with the box cone's t weighted like the other rows of the cone), a
 * working warm start gives at most 1.6e-2: its first iterate moves off the
 * solution by up to 1.6% of a cold start's distance where the solve it starts
 * from stopped mid-swing (see hs21_tiny_qp). With BOX_T_RY_FACTOR at 1e4,
 * its value in double precision (Spark gcc, double and single precision,
 * QDLDL and CG, the same 2,000 scales), it gives at most 9.6e-3.
 * Deliberately broken warm starts in SCS give at least 1 when it is ignored
 * (a cold start), 0.95 with x left at zero, 1.1 with s dropped from
 * v = (x, y + R^{-1} s, 1), 4.0 when v is not normalized, 750 with s scaled
 * by R instead of R^{-1}, and 0.12 and 0.24 with y dropped or negated, the
 * subtlest since hs21's y is small. A block where a cold start's first
 * iterate is already at sol cannot tell the two apart; it fails the check
 * unless the warm start's is there too. */
#define WARM_START_FIRST_ITER_TOL (5e-2)

/* max(a, b), except that it is NaN if a or b is NaN, so once a running
 * maximum has met a NaN it stays NaN whatever is compared after it */
static scs_float warm_start_nan_max(scs_float a, scs_float b) {
  if (isnan(a)) {
    return a;
  }
  return (isnan(b) || b > a) ? b : a;
}

scs_float warm_start_first_iterate_ratio(const ScsData *d, const ScsCone *k,
                                         const ScsSettings *stgs,
                                         const ScsSolution *sol) {
  scs_int n = d->n, m = d->m, i, b, warm, len[3];
  scs_float dist[2][3] = {{0.}}, ratio = 0., e;
  const scs_float *ref[3];
  scs_float *it[3];
  ScsSettings st = *stgs;
  ScsSolution first = {SCS_NULL, SCS_NULL, SCS_NULL};
  ScsInfo info = {0};
  first.x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  first.y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  first.s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  ref[0] = sol->x;
  ref[1] = sol->y;
  ref[2] = sol->s;
  it[0] = first.x;
  it[1] = first.y;
  it[2] = first.s;
  len[0] = n;
  len[1] = len[2] = m;
  st.max_iters = 1;
  for (warm = 0; warm < 2; ++warm) {
    st.warm_start = warm;
    for (b = 0; b < 3; ++b) {
      memcpy(it[b], ref[b], len[b] * sizeof(scs_float));
    }
    scs(d, k, &st, &first, &info);
    /* a NaN anywhere in a block makes its distance NaN */
    for (b = 0; b < 3; ++b) {
      for (i = 0; i < len[b]; ++i) {
        dist[warm][b] =
            warm_start_nan_max(dist[warm][b], ABS(it[b][i] - ref[b][i]));
      }
    }
  }
  /* a NaN distance in any block, warm or cold, makes the ratio NaN, which
   * fails a ratio <= tolerance check */
  for (b = 0; b < 3; ++b) {
    e = (dist[1][b] == 0. && !isnan(dist[0][b]))
            ? 0.
            : dist[1][b] / dist[0][b];
    ratio = warm_start_nan_max(ratio, e);
  }
  scs_free(first.x);
  scs_free(first.y);
  scs_free(first.s);
  return ratio;
}

#endif
