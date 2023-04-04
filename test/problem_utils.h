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

#define _MAX_RAND_VAL (1073741823) /* 2^30 - 1 */

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
  ran_start(seed);
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
    mu_assert_less("Primal residual ERROR", ABS(res_pri - info->res_pri),
                   1e-10);
    mu_assert_less("Dual residual ERROR", ABS(res_dual - info->res_dual),
                   1e-10);
    mu_assert_less("Gap ERROR", ABS(gap - info->gap), 1e-7 * (1 + ABS(gap)));
    mu_assert_less("Primal obj ERROR", ABS(pobj - info->pobj),
                   1e-9 * (1 + ABS(pobj)));
    mu_assert_less("Dual obj ERROR", ABS(dobj - info->dobj),
                   1e-9 * (1 + ABS(dobj)));
    /* slightly looser tol */
    mu_assert_less("Complementary slackness ERROR", ABS(sty),
                   5e-8 * MAX(NORM(s, m), NORM(y, m)));
    mu_assert_less("s cone dist ERROR", ABS(sdist), 1e-5);
    mu_assert_less("y cone dist ERROR", ABS(ydist), 1e-5);

    mu_assert_less("Primal feas ERROR", res_pri,
                   stgs->eps_abs + stgs->eps_rel * prl);
    mu_assert_less("Dual feas ERROR", res_dual,
                   stgs->eps_abs + stgs->eps_rel * drl);
    mu_assert_less("Gap feas ERROR", gap, stgs->eps_abs + stgs->eps_rel * grl);

  } else if (status == SCS_INFEASIBLE) {
    mu_assert_less("Infeas ERROR", ABS(res_infeas - info->res_infeas), 1e-8);
    mu_assert_less("bty ERROR", ABS(bty + 1), 1e-12);
    mu_assert_less("y cone dist ERROR", ABS(ydist), 1e-5);
    mu_assert_less("Infeas invalid ERROR", res_infeas, stgs->eps_infeas);

  } else if (status == SCS_UNBOUNDED) {
    mu_assert_less("Unbdd_a ERROR", ABS(res_unbdd_a - info->res_unbdd_a), 1e-8);
    mu_assert_less("Unbdd_p ERROR", ABS(res_unbdd_p - info->res_unbdd_p), 1e-8);
    mu_assert_less("ctx ERROR", ABS(ctx + 1), 1e-12);
    mu_assert_less("s cone dist ERROR", ABS(sdist), 1e-5);
    mu_assert_less("Unbounded P invalid ERROR", res_unbdd_p, stgs->eps_infeas);
    mu_assert_less("Unbounded A invalid ERROR", res_unbdd_a, stgs->eps_infeas);

  } else {
    return "INVALID STATUS";
  }
  return 0;
}

#endif
