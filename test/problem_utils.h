#ifndef PUTILS_H_GUARD
#define PUTILS_H_GUARD

#include "amatrix.h"
#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "util.h"
#include "rng.h"

#define PI (3.141592654)
#ifdef DLONG
#ifdef _WIN64
/* this is a Microsoft extension, but also works with min_g_w-w64 */
#define INTRW "%I64d"
#else
#define INTRW "%ld"
#endif
#else
#define INTRW "%i"
#endif

/* uniform random number in [-1,1] */
static scs_float rand_scs_float(void) {
  return 2 * (((scs_float)ran_arr_next()) / RAND_MAX) - 1;
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
  tmp_cone_work = SCS(init_cone)(k, SCS_NULL, m);
  SCS(proj_dual_cone)(y, k, tmp_cone_work, 0);
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


const char * verify_solution_correct(ScsData * d, ScsCone * k, ScsInfo *info, ScsSolution *sol) {
  scs_int n = d->n, m = d->m;
  scs_float *x = sol->x;
  scs_float *y = sol->y;
  scs_float *s = sol->s;

  scs_float *c = d->c;
  scs_float *b = d->b;

  scs_float *primal = scs_calloc(m, sizeof(scs_float));
  scs_float *dual = scs_calloc(n, sizeof(scs_float));
  scs_float *px = scs_calloc(n, sizeof(scs_float));
  
  scs_float res_pri, res_dual, res_infeas, res_unbdd_a, res_unbdd_p;
  scs_float ctx, bty, xt_p_x, gap, pobj, dobj;

  /**************** PRIMAL *********************/
  memset(primal, 0, m * sizeof(scs_float));
  /* Ax */
  SCS(accum_by_a)(d->A, SCS_NULL, x, primal);
  /* Ax + s */
  SCS(add_scaled_array)(primal, s, m, 1.);

  res_infeas = NORM(primal, m);

  /* Ax + s - b * tau */
  SCS(add_scaled_array)(primal, b, m, -1.0);

  res_pri = NORM(primal, m); 

  /**************** DUAL *********************/
  memset(px, 0, n * sizeof(scs_float));
  if (d->P) {
    /* px = Px */
    SCS(accum_by_p)(d->P, SCS_NULL, x, px);
    xt_p_x = SCS(dot)(px, x, n);
    res_unbdd_p = NORM(px, n); 
  } else{
    xt_p_x = 0;
    res_unbdd_p = 0;
  }

  memset(dual, 0, n * sizeof(scs_float));
  /* aty = A'y */
  SCS(accum_by_atrans)(d->A, SCS_NULL, y, dual);
  res_unbdd_a = NORM(dual, n);  

  /* Px + A'y */
  SCS(add_scaled_array)(dual, px, n, 1.);
  /* Px + A'y + c */
  SCS(add_scaled_array)(dual, c, n, 1.0);
  
  res_dual = NORM(dual, n);

  /**************** OTHERS *****************/
  bty = SCS(dot)(y, b, m);
  ctx = SCS(dot)(x, c, n);

  gap = ABS(xt_p_x + ctx + bty);
  pobj = xt_p_x / 2. + ctx;
  dobj = -xt_p_x / 2. - bty;

  /**************** ASSERTS *****************/

  mu_assert("Primal residual wrong", ABS(res_pri - info->res_pri) < 1e-10);
  mu_assert("Dual residual wrong", ABS(res_dual - info->res_dual) < 1e-10);
  mu_assert("Infeas wrong", ABS(res_infeas - info->res_infeas) < 1e-10);
  mu_assert("Unbdd_a wrong", ABS(res_unbdd_a - info->res_unbdd_a) < 1e-10);
  mu_assert("Unbdd_p wrong", ABS(res_unbdd_p - info->res_unbdd_p) < 1e-10);
  mu_assert("Gap wrong", ABS(gap - info->gap) < 1e-10);
  mu_assert("Primal obj wrong", ABS(pobj - info->pobj) < 1e-10);
  mu_assert("Dual obj wrong", ABS(dobj - info->dobj) < 1e-10);

  /**************** CLEANUP *****************/
  scs_free(primal);
  scs_free(dual);
  scs_free(px);
  return 0;
}


#endif
