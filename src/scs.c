#include "scs.h"

#include "aa.h"
#include "ctrlc.h"
#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "normalize.h"
#include "rw.h"
#include "util.h"

SCS(timer) global_timer;

/* printing header */
static const char *HEADER[] = {
    " Iter ",    " pri res ", " dua res ", " rel gap ",
    " pri obj ", " dua obj ", " kap/tau ", " time (s)",
};
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 8;
static const scs_int LINE_LEN = 76;

static scs_int scs_isnan(scs_float x) { return (x == NAN || x != x); }

static void free_work(ScsWork *w) {
  if (w) {
    scs_free(w->u);
    scs_free(w->u_best);
    scs_free(w->u_t);
    scs_free(w->u_prev);
    /* Don't need these because u*, v* are contiguous in mem
      scs_free(w->v);
      scs_free(w->v_best);
      scs_free(w->v_prev);
    */
    scs_free(w->h);
    scs_free(w->g);
    scs_free(w->b);
    scs_free(w->c);
    scs_free(w->pr);
    scs_free(w->dr);
    if (w->scal) {
      scs_free(w->scal->D);
      scs_free(w->scal->E);
      scs_free(w->scal);
    }
    scs_free(w);
  }
}

static void print_init_header(const ScsData *d, const ScsCone *k) {
  scs_int i;
  ScsSettings *stgs = d->stgs;
  char *cone_str = SCS(get_cone_header)(k);
  char *lin_sys_method = SCS(get_lin_sys_method)(d->A, d->stgs);
#ifdef USE_LAPACK
  scs_int acceleration_lookback = stgs->acceleration_lookback;
#else
  scs_int acceleration_lookback = 0;
#endif
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf(
      "\n\tSCS v%s - Splitting Conic Solver\n\t(c) Brendan "
      "O'Donoghue, Stanford University, 2012\n",
      SCS(version)());
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  if (lin_sys_method) {
    scs_printf("Lin-sys: %s\n", lin_sys_method);
    scs_free(lin_sys_method);
  }
  if (stgs->normalize) {
    scs_printf(
        "eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i, "
        "scale = %2.2f\nacceleration_lookback = %i, rho_x = %.2e\n",
        stgs->eps, stgs->alpha, (int)stgs->max_iters, (int)stgs->normalize,
        stgs->scale, (int)acceleration_lookback, stgs->rho_x);
  } else {
    scs_printf(
        "eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i\n"
        "acceleration_lookback = %i, rho_x = %.2e\n",
        stgs->eps, stgs->alpha, (int)stgs->max_iters, (int)stgs->normalize,
        (int)acceleration_lookback, stgs->rho_x);
  }
  scs_printf("Variables n = %i, constraints m = %i\n", (int)d->n, (int)d->m);
  scs_printf("%s", cone_str);
  scs_free(cone_str);
#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

static void populate_on_failure(scs_int m, scs_int n, ScsSolution *sol,
                                ScsInfo *info, scs_int status_val,
                                const char *msg) {
  if (info) {
    info->rel_gap = NAN;
    info->res_pri = NAN;
    info->res_dual = NAN;
    info->pobj = NAN;
    info->dobj = NAN;
    info->iter = -1;
    info->status_val = status_val;
    info->solve_time = NAN;
    strcpy(info->status, msg);
  }
  if (sol) {
    if (n > 0) {
      if (!sol->x) {
        sol->x = (scs_float *)scs_malloc(sizeof(scs_float) * n);
      }
      SCS(scale_array)(sol->x, NAN, n);
    }
    if (m > 0) {
      if (!sol->y) {
        sol->y = (scs_float *)scs_malloc(sizeof(scs_float) * m);
      }
      SCS(scale_array)(sol->y, NAN, m);
      if (!sol->s) {
        sol->s = (scs_float *)scs_malloc(sizeof(scs_float) * m);
      }
      SCS(scale_array)(sol->s, NAN, m);
    }
  }
}

static scs_int failure(ScsWork *w, scs_int m, scs_int n, ScsSolution *sol,
                       ScsInfo *info, scs_int stint, const char *msg,
                       const char *ststr) {
  scs_int status = stint;
  populate_on_failure(m, n, sol, info, status, ststr);
  scs_printf("Failure:%s\n", msg);
  scs_end_interrupt_listener();
  return status;
}

static void warm_start_vars(ScsWork *w, const ScsSolution *sol) {
  scs_int i, n = w->n, m = w->m;
  memset(w->v, 0, n * sizeof(scs_float));
  memcpy(w->u, sol->x, n * sizeof(scs_float));
  memcpy(&(w->u[n]), sol->y, m * sizeof(scs_float));
  memcpy(&(w->v[n]), sol->s, m * sizeof(scs_float));
  w->u[n + m] = 1.0;
  w->v[n + m] = 0.0;
#ifndef NOVALIDATE
  for (i = 0; i < n + m + 1; ++i) {
    if (scs_isnan(w->u[i])) {
      w->u[i] = 0;
    }
    if (scs_isnan(w->v[i])) {
      w->v[i] = 0;
    }
  }
#endif
  if (w->stgs->normalize) {
    SCS(normalize_warm_start)(w);
  }
}

static scs_float calc_primal_resid(ScsWork *w, const scs_float *x,
                                   const scs_float *s, const scs_float tau,
                                   scs_float *nm_axs) {
  scs_int i;
  scs_float pres = 0, scale, *pr = w->pr;
  *nm_axs = 0;
  memset(pr, 0, w->m * sizeof(scs_float));
  SCS(accum_by_a)(w->A, w->p, x, pr);
  SCS(add_scaled_array)(pr, s, w->m, 1.0); /* pr = Ax + s */
  for (i = 0; i < w->m; ++i) {
    scale = w->stgs->normalize ? w->scal->D[i] / (w->sc_b * w->stgs->scale) : 1;
    scale = scale * scale;
    *nm_axs += (pr[i] * pr[i]) * scale;
    pres += (pr[i] - w->b[i] * tau) * (pr[i] - w->b[i] * tau) * scale;
  }
  *nm_axs = SQRTF(*nm_axs);
  return SQRTF(pres); /* SCS(norm)(Ax + s - b * tau) */
}

static scs_float calc_dual_resid(ScsWork *w, const scs_float *y,
                                 const scs_float tau, scs_float *nm_a_ty) {
  scs_int i;
  scs_float dres = 0, scale, *dr = w->dr;
  *nm_a_ty = 0;
  memset(dr, 0, w->n * sizeof(scs_float));
  SCS(accum_by_atrans)(w->A, w->p, y, dr); /* dr = A'y */
  for (i = 0; i < w->n; ++i) {
    scale = w->stgs->normalize ? w->scal->E[i] / (w->sc_c * w->stgs->scale) : 1;
    scale = scale * scale;
    *nm_a_ty += (dr[i] * dr[i]) * scale;
    dres += (dr[i] + w->c[i] * tau) * (dr[i] + w->c[i] * tau) * scale;
  }
  *nm_a_ty = SQRTF(*nm_a_ty);
  return SQRTF(dres); /* SCS(norm)(A'y + c * tau) */
}

/* calculates un-normalized quantities */
static void calc_residuals(ScsWork *w, ScsResiduals *r, scs_int iter) {
  scs_float *x = w->u, *y = &(w->u[w->n]), *s = &(w->v[w->n]);
  scs_float nmpr_tau, nmdr_tau, nm_axs_tau, nm_a_ty_tau, ct_x, bt_y;
  scs_int n = w->n, m = w->m;

  /* checks if the residuals are unchanged by checking iteration */
  if (r->last_iter == iter) {
    return;
  }
  r->last_iter = iter;

  r->tau = ABS(w->u[n + m]);
  r->kap = ABS(w->v[n + m]) /
           (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);

  nmpr_tau = calc_primal_resid(w, x, s, r->tau, &nm_axs_tau);
  nmdr_tau = calc_dual_resid(w, y, r->tau, &nm_a_ty_tau);

  r->bt_y_by_tau =
      SCS(dot)(y, w->b, m) /
      (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);
  r->ct_x_by_tau =
      SCS(dot)(x, w->c, n) /
      (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);

  r->res_infeas =
      r->bt_y_by_tau < 0 ? w->nm_b * nm_a_ty_tau / -r->bt_y_by_tau : NAN;
  r->res_unbdd =
      r->ct_x_by_tau < 0 ? w->nm_c * nm_axs_tau / -r->ct_x_by_tau : NAN;

  bt_y = SAFEDIV_POS(r->bt_y_by_tau, r->tau);
  ct_x = SAFEDIV_POS(r->ct_x_by_tau, r->tau);

  r->res_pri = SAFEDIV_POS(nmpr_tau / (1 + w->nm_b), r->tau);
  r->res_dual = SAFEDIV_POS(nmdr_tau / (1 + w->nm_c), r->tau);
  r->rel_gap = ABS(ct_x + bt_y) / (1 + ABS(ct_x) + ABS(bt_y));
}

static void cold_start_vars(ScsWork *w) {
  scs_int l = w->n + w->m + 1;
  memset(w->u, 0, l * sizeof(scs_float));
  memset(w->v, 0, l * sizeof(scs_float));
  w->u[l - 1] = SQRTF((scs_float)l);
  w->v[l - 1] = SQRTF((scs_float)l);
}

/* status < 0 indicates failure */
static scs_int project_lin_sys(ScsWork *w, scs_int iter) {
  /* ut = u + v */

  scs_int n = w->n, m = w->m, l = n + m + 1, status;
  memcpy(w->u_t, w->u, l * sizeof(scs_float));
  SCS(add_scaled_array)(w->u_t, w->v, l, 1.0);

  SCS(scale_array)(w->u_t, w->stgs->rho_x, n);

  SCS(add_scaled_array)(w->u_t, w->h, l - 1, -w->u_t[l - 1]);
  SCS(add_scaled_array)
  (w->u_t, w->h, l - 1, -SCS(dot)(w->u_t, w->g, l - 1) / (w->g_th + 1));
  SCS(scale_array)(&(w->u_t[n]), -1, m);

  status = SCS(solve_lin_sys)(w->A, w->stgs, w->p, w->u_t, w->u, iter);

  w->u_t[l - 1] += SCS(dot)(w->u_t, w->h, l - 1);

  return status;
}

static void update_dual_vars(ScsWork *w) {
  scs_int i, n = w->n, l = n + w->m + 1;
  /* this does not relax 'x' variable */
  for (i = n; i < l; ++i) {
    w->v[i] += (w->u[i] - w->stgs->alpha * w->u_t[i] -
                (1.0 - w->stgs->alpha) * w->u_prev[i]);
  }
}

/* status < 0 indicates failure */
static scs_int project_cones(ScsWork *w, const ScsCone *k, scs_int iter) {
  scs_int i, n = w->n, l = n + w->m + 1, status;
  /* this does not relax 'x' variable */
  for (i = 0; i < n; ++i) {
    w->u[i] = w->u_t[i] - w->v[i];
  }
  for (i = n; i < l; ++i) {
    w->u[i] = w->stgs->alpha * w->u_t[i] + (1 - w->stgs->alpha) * w->u_prev[i] -
              w->v[i];
  }
  /* u = [x;y;tau] */
  status =
      SCS(proj_dual_cone)(&(w->u[n]), k, w->cone_work, &(w->u_prev[n]), iter);
  if (w->u[l - 1] < 0.0) {
    w->u[l - 1] = 0.0;
  }

  return status;
}

static scs_int indeterminate(ScsWork *w, ScsSolution *sol, ScsInfo *info) {
  strcpy(info->status, "Indeterminate");
  SCS(scale_array)(sol->x, NAN, w->n);
  SCS(scale_array)(sol->y, NAN, w->m);
  SCS(scale_array)(sol->s, NAN, w->m);
  return SCS_INDETERMINATE;
}

static void sety(ScsWork *w, ScsSolution *sol) {
  if (!sol->y) {
    sol->y = (scs_float *)scs_malloc(sizeof(scs_float) * w->m);
  }
  memcpy(sol->y, &(w->u[w->n]), w->m * sizeof(scs_float));
}

static void sets(ScsWork *w, ScsSolution *sol) {
  if (!sol->s) {
    sol->s = (scs_float *)scs_malloc(sizeof(scs_float) * w->m);
  }
  memcpy(sol->s, &(w->v[w->n]), w->m * sizeof(scs_float));
}

static void setx(ScsWork *w, ScsSolution *sol) {
  if (!sol->x) {
    sol->x = (scs_float *)scs_malloc(sizeof(scs_float) * w->n);
  }
  memcpy(sol->x, w->u, w->n * sizeof(scs_float));
}

static scs_float get_max_residual(ScsResiduals *r) {
  return MAX(r->rel_gap, MAX(r->res_pri, r->res_dual));
}

static void copy_from_best_iterate(ScsWork *w) {
  memcpy(w->u, w->u_best, (w->m + w->n + 1) * sizeof(scs_float));
  memcpy(w->v, w->v_best, (w->m + w->n + 1) * sizeof(scs_float));
}

static scs_int solved(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                      ScsResiduals *r, scs_int iter) {
  if (w->best_max_residual < get_max_residual(r)) {
    r->last_iter = -1; /* Forces residual recomputation. */
    copy_from_best_iterate(w);
    calc_residuals(w, r, iter);
    setx(w, sol);
    sety(w, sol);
    sets(w, sol);
  }
  SCS(scale_array)(sol->x, SAFEDIV_POS(1.0, r->tau), w->n);
  SCS(scale_array)(sol->y, SAFEDIV_POS(1.0, r->tau), w->m);
  SCS(scale_array)(sol->s, SAFEDIV_POS(1.0, r->tau), w->m);
  if (info->status_val == 0) {
    strcpy(info->status, "Solved/Inaccurate");
    return SCS_SOLVED_INACCURATE;
  }
  strcpy(info->status, "Solved");
  return SCS_SOLVED;
}

static scs_int infeasible(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                          scs_float bt_y) {
  SCS(scale_array)(sol->y, -1 / bt_y, w->m);
  SCS(scale_array)(sol->x, NAN, w->n);
  SCS(scale_array)(sol->s, NAN, w->m);
  if (info->status_val == 0) {
    strcpy(info->status, "Infeasible/Inaccurate");
    return SCS_INFEASIBLE_INACCURATE;
  }
  strcpy(info->status, "Infeasible");
  return SCS_INFEASIBLE;
}

static scs_int unbounded(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                         scs_float ct_x) {
  SCS(scale_array)(sol->x, -1 / ct_x, w->n);
  SCS(scale_array)(sol->s, -1 / ct_x, w->m);
  SCS(scale_array)(sol->y, NAN, w->m);
  if (info->status_val == 0) {
    strcpy(info->status, "Unbounded/Inaccurate");
    return SCS_UNBOUNDED_INACCURATE;
  }
  strcpy(info->status, "Unbounded");
  return SCS_UNBOUNDED;
}

static scs_int is_solved_status(scs_int status) {
  return status == SCS_SOLVED || status == SCS_SOLVED_INACCURATE;
}

static scs_int is_infeasible_status(scs_int status) {
  return status == SCS_INFEASIBLE || status == SCS_INFEASIBLE_INACCURATE;
}

static scs_int is_unbounded_status(scs_int status) {
  return status == SCS_UNBOUNDED || status == SCS_UNBOUNDED_INACCURATE;
}

static void get_info(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                     ScsResiduals *r, scs_int iter) {
  info->iter = iter;
  info->res_infeas = r->res_infeas;
  info->res_unbdd = r->res_unbdd;
  if (is_solved_status(info->status_val)) {
    info->rel_gap = r->rel_gap;
    info->res_pri = r->res_pri;
    info->res_dual = r->res_dual;
    info->pobj = r->ct_x_by_tau / r->tau;
    info->dobj = -r->bt_y_by_tau / r->tau;
  } else if (is_unbounded_status(info->status_val)) {
    info->rel_gap = NAN;
    info->res_pri = NAN;
    info->res_dual = NAN;
    info->pobj = -INFINITY;
    info->dobj = -INFINITY;
  } else if (is_infeasible_status(info->status_val)) {
    info->rel_gap = NAN;
    info->res_pri = NAN;
    info->res_dual = NAN;
    info->pobj = INFINITY;
    info->dobj = INFINITY;
  }
}

/* sets solutions, re-scales by inner prods if infeasible or unbounded */
static void get_solution(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                         ScsResiduals *r, scs_int iter) {
  scs_int l = w->n + w->m + 1;
  calc_residuals(w, r, iter);
  setx(w, sol);
  sety(w, sol);
  sets(w, sol);
  if (info->status_val == SCS_UNFINISHED) {
    /* not yet converged, take best guess */
    if (r->tau > INDETERMINATE_TOL && r->tau > r->kap) {
      info->status_val = solved(w, sol, info, r, iter);
    } else if (SCS(norm)(w->u, l) < INDETERMINATE_TOL * SQRTF((scs_float)l)) {
      info->status_val = indeterminate(w, sol, info);
    } else if (r->bt_y_by_tau < r->ct_x_by_tau) {
      info->status_val = infeasible(w, sol, info, r->bt_y_by_tau);
    } else {
      info->status_val = unbounded(w, sol, info, r->ct_x_by_tau);
    }
  } else if (is_solved_status(info->status_val)) {
    info->status_val = solved(w, sol, info, r, iter);
  } else if (is_infeasible_status(info->status_val)) {
    info->status_val = infeasible(w, sol, info, r->bt_y_by_tau);
  } else {
    info->status_val = unbounded(w, sol, info, r->ct_x_by_tau);
  }
  if (w->stgs->normalize) {
    SCS(un_normalize_sol)(w, sol);
  }
  get_info(w, sol, info, r, iter);
}

static void print_summary(ScsWork *w, scs_int i, ScsResiduals *r,
                          SCS(timer) * solve_timer) {
  scs_printf("%*i|", (int)strlen(HEADER[0]), (int)i);
  scs_printf("%*.2e ", (int)HSPACE, r->res_pri);
  scs_printf("%*.2e ", (int)HSPACE, r->res_dual);
  scs_printf("%*.2e ", (int)HSPACE, r->rel_gap);
  scs_printf("%*.2e ", (int)HSPACE, SAFEDIV_POS(r->ct_x_by_tau, r->tau));
  scs_printf("%*.2e ", (int)HSPACE, SAFEDIV_POS(-r->bt_y_by_tau, r->tau));
  scs_printf("%*.2e ", (int)HSPACE, SAFEDIV_POS(r->kap, r->tau));
  scs_printf("%*.2e ", (int)HSPACE, SCS(tocq)(solve_timer) / 1e3);
  scs_printf("\n");

#if EXTRA_VERBOSE > 0
  scs_printf("Norm u = %4f, ", SCS(norm)(w->u, w->n + w->m + 1));
  scs_printf("Norm u_t = %4f, ", SCS(norm)(w->u_t, w->n + w->m + 1));
  scs_printf("Norm v = %4f, ", SCS(norm)(w->v, w->n + w->m + 1));
  scs_printf("tau = %4f, ", w->u[w->n + w->m]);
  scs_printf("kappa = %4f, ", w->v[w->n + w->m]);
  scs_printf("|u - u_prev| = %1.2e, ",
             SCS(norm_diff)(w->u, w->u_prev, w->n + w->m + 1));
  scs_printf("|u - u_t| = %1.2e, ",
             SCS(norm_diff)(w->u, w->u_t, w->n + w->m + 1));
  scs_printf("res_infeas = %1.2e, ", r->res_infeas);
  scs_printf("res_unbdd = %1.2e\n", r->res_unbdd);
#endif

#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

static void print_header(ScsWork *w, const ScsCone *k) {
  scs_int i;
  if (w->stgs->warm_start) {
    scs_printf("SCS using variable warm-starting\n");
  }
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  for (i = 0; i < HEADER_LEN - 1; ++i) {
    scs_printf("%s|", HEADER[i]);
  }
  scs_printf("%s\n", HEADER[HEADER_LEN - 1]);
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

static scs_float get_dual_cone_dist(const scs_float *y, const ScsCone *k,
                                    ScsConeWork *c, scs_int m) {
  scs_float dist;
  scs_float *t = (scs_float *)scs_malloc(sizeof(scs_float) * m);
  memcpy(t, y, m * sizeof(scs_float));
  SCS(proj_dual_cone)(t, k, c, SCS_NULL, -1);
  dist = SCS(norm_inf_diff)(t, y, m);
#if EXTRA_VERBOSE > 0
  SCS(print_array)(y, m, "y");
  SCS(print_array)(t, m, "proj_y");
  scs_printf("dist = %4f\n", dist);
#endif
  scs_free(t);
  return dist;
}

/* via moreau */
static scs_float get_pri_cone_dist(const scs_float *s, const ScsCone *k,
                                   ScsConeWork *c, scs_int m) {
  scs_float dist;
  scs_float *t = (scs_float *)scs_malloc(sizeof(scs_float) * m);
  memcpy(t, s, m * sizeof(scs_float));
  SCS(scale_array)(t, -1.0, m);
  SCS(proj_dual_cone)(t, k, c, SCS_NULL, -1);
  dist = SCS(norm_inf)(t, m); /* ||s - Pi_c(s)|| = ||Pi_c*(-s)|| */
#if EXTRA_VERBOSE > 0
  SCS(print_array)(s, m, "s");
  SCS(print_array)(t, m, "(s - proj_s)");
  scs_printf("dist = %4f\n", dist);
#endif
  scs_free(t);
  return dist;
}

static char *get_accel_summary(ScsInfo *info, scs_float total_accel_time) {
  char *str = (char *)scs_malloc(sizeof(char) * 64);
  sprintf(str, "\tAcceleration: avg step time: %1.2es\n",
          total_accel_time / (info->iter + 1) / 1e3);
  return str;
}

static void print_footer(const ScsData *d, const ScsCone *k, ScsSolution *sol,
                         ScsWork *w, ScsInfo *info,
                         scs_float total_accel_time) {
  scs_int i;
  char *lin_sys_str = SCS(get_lin_sys_summary)(w->p, info);
  char *cone_str = SCS(get_cone_summary)(info, w->cone_work);
  char *accel_str = get_accel_summary(info, total_accel_time);
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\nStatus: %s\n", info->status);
  if (info->iter == w->stgs->max_iters) {
    scs_printf(
        "Hit max_iters, solution may be inaccurate, returning best found "
        "solution.\n");
  }
  scs_printf("Timing: Solve time: %1.2es\n", info->solve_time / 1e3);

  if (lin_sys_str) {
    scs_printf("%s", lin_sys_str);
    scs_free(lin_sys_str);
  }

  if (cone_str) {
    scs_printf("%s", cone_str);
    scs_free(cone_str);
  }

  if (accel_str) {
    scs_printf("%s", accel_str);
    scs_free(accel_str);
  }

  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");

  if (is_infeasible_status(info->status_val)) {
    scs_printf("Certificate of primal infeasibility:\n");
    scs_printf("dist(y, K*) = %.4e\n",
               get_dual_cone_dist(sol->y, k, w->cone_work, d->m));
    scs_printf("|A'y|_2 * |b|_2 = %.4e\n", info->res_infeas);
    scs_printf("b'y = %.4f\n", SCS(dot)(d->b, sol->y, d->m));
  } else if (is_unbounded_status(info->status_val)) {
    scs_printf("Certificate of dual infeasibility:\n");
    scs_printf("dist(s, K) = %.4e\n",
               get_pri_cone_dist(sol->s, k, w->cone_work, d->m));
    scs_printf("|Ax + s|_2 * |c|_2 = %.4e\n", info->res_unbdd);
    scs_printf("c'x = %.4f\n", SCS(dot)(d->c, sol->x, d->n));
  } else {
    scs_printf("Error metrics:\n");
    scs_printf("dist(s, K) = %.4e, dist(y, K*) = %.4e, s'y/|s||y| = %.4e\n",
               get_pri_cone_dist(sol->s, k, w->cone_work, d->m),
               get_dual_cone_dist(sol->y, k, w->cone_work, d->m),
               SCS(dot)(sol->s, sol->y, d->m) / SCS(norm)(sol->s, d->m) /
                   SCS(norm)(sol->y, d->m));
    scs_printf("primal res: |Ax + s - b|_2 / (1 + |b|_2) = %.4e\n",
               info->res_pri);
    scs_printf("dual res:   |A'y + c|_2 / (1 + |c|_2) = %.4e\n",
               info->res_dual);
    scs_printf("rel gap:    |c'x + b'y| / (1 + |c'x| + |b'y|) = %.4e\n",
               info->rel_gap);
    for (i = 0; i < LINE_LEN; ++i) {
      scs_printf("-");
    }
    scs_printf("\n");
    scs_printf("c'x = %.4f, -b'y = %.4f\n", info->pobj, info->dobj);
  }
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("=");
  }
  scs_printf("\n");
#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

static scs_int has_converged(ScsWork *w, ScsResiduals *r, scs_int iter) {
  scs_float eps = w->stgs->eps;
  if (isless(r->res_pri, eps) && isless(r->res_dual, eps) &&
      isless(r->rel_gap, eps)) {
    return SCS_SOLVED;
  }
  /* Add iter > 0 to avoid strange edge case where infeasible point found
   * right at start of run `out/demo_SOCP_indirect 2 0.1 0.3 1506264403` */
  if (isless(r->res_unbdd, eps) && iter > 0) {
    return SCS_UNBOUNDED;
  }
  if (isless(r->res_infeas, eps) && iter > 0) {
    return SCS_INFEASIBLE;
  }
  return 0;
}

static scs_int validate(const ScsData *d, const ScsCone *k) {
  ScsSettings *stgs = d->stgs;
  if (d->m <= 0 || d->n <= 0) {
    scs_printf("m and n must both be greater than 0; m = %li, n = %li\n",
               (long)d->m, (long)d->n);
    return -1;
  }
  if (d->m < d->n) {
    scs_printf("WARN: m less than n, problem likely degenerate\n");
    /* return -1; */
  }
  if (SCS(validate_lin_sys)(d->A) < 0) {
    scs_printf("invalid linear system input data\n");
    return -1;
  }
  if (SCS(validate_cones)(d, k) < 0) {
    scs_printf("cone validation error\n");
    return -1;
  }
  if (stgs->max_iters <= 0) {
    scs_printf("max_iters must be positive\n");
    return -1;
  }
  if (stgs->eps <= 0) {
    scs_printf("eps tolerance must be positive\n");
    return -1;
  }
  if (stgs->alpha <= 0 || stgs->alpha >= 2) {
    scs_printf("alpha must be in (0,2)\n");
    return -1;
  }
  if (stgs->rho_x <= 0) {
    scs_printf("rho_x must be positive (1e-3 works well).\n");
    return -1;
  }
  if (stgs->scale <= 0) {
    scs_printf("scale must be positive (1 works well).\n");
    return -1;
  }
  return 0;
}

static ScsWork *init_work(const ScsData *d, const ScsCone *k) {
  ScsWork *w = (ScsWork *)scs_calloc(1, sizeof(ScsWork));
  scs_int l = d->n + d->m + 1;
  if (d->stgs->verbose) {
    print_init_header(d, k);
  }
  if (!w) {
    scs_printf("ERROR: allocating work failure\n");
    return SCS_NULL;
  }
  /* get settings and dims from data struct */
  w->stgs = d->stgs;
  w->m = d->m;
  w->n = d->n;
  w->best_max_residual = INFINITY;
  /* allocate workspace: */
  /* u* include v* values */
  w->u = (scs_float *)scs_malloc(2 * l * sizeof(scs_float));
  w->u_best = (scs_float *)scs_malloc(2 * l * sizeof(scs_float));
  w->u_t = (scs_float *)scs_malloc(l * sizeof(scs_float));
  w->u_prev = (scs_float *)scs_malloc(2 * l * sizeof(scs_float));
  w->h = (scs_float *)scs_malloc((l - 1) * sizeof(scs_float));
  w->g = (scs_float *)scs_malloc((l - 1) * sizeof(scs_float));
  w->pr = (scs_float *)scs_malloc(d->m * sizeof(scs_float));
  w->dr = (scs_float *)scs_malloc(d->n * sizeof(scs_float));
  w->b = (scs_float *)scs_malloc(d->m * sizeof(scs_float));
  w->c = (scs_float *)scs_malloc(d->n * sizeof(scs_float));
  if (!w->u || !w->u_t || !w->u_prev || !w->h || !w->g || !w->pr || !w->dr ||
      !w->b || !w->c) {
    scs_printf("ERROR: work memory allocation failure\n");
    return SCS_NULL;
  }
  /* make u,v and u_prev,v_prev contiguous in memory */
  w->v = &(w->u[l]);
  w->v_best = &(w->u_best[l]);
  w->v_prev = &(w->u_prev[l]);
  w->A = d->A;
  if (w->stgs->normalize) {
#ifdef COPYAMATRIX
    if (!SCS(copy_a_matrix)(&(w->A), d->A)) {
      scs_printf("ERROR: copy A matrix failed\n");
      return SCS_NULL;
    }
#endif
    w->scal = (ScsScaling *)scs_malloc(sizeof(ScsScaling));
    SCS(normalize_a)(w->A, w->stgs, k, w->scal);
#if EXTRA_VERBOSE > 0
    SCS(print_array)(w->scal->D, d->m, "D");
    scs_printf("SCS(norm) D = %4f\n", SCS(norm)(w->scal->D, d->m));
    SCS(print_array)(w->scal->E, d->n, "E");
    scs_printf("SCS(norm) E = %4f\n", SCS(norm)(w->scal->E, d->n));
#endif
  } else {
    w->scal = SCS_NULL;
  }
  if (!(w->cone_work = SCS(init_cone)(k))) {
    scs_printf("ERROR: init_cone failure\n");
    return SCS_NULL;
  }
  if (!(w->p = SCS(init_lin_sys_work)(w->A, w->stgs))) {
    scs_printf("ERROR: init_lin_sys_work failure\n");
    return SCS_NULL;
  }
  if (!(w->accel =
            aa_init(2 * (w->m + w->n + 1), ABS(w->stgs->acceleration_lookback),
                    w->stgs->acceleration_lookback >= 0))) {
    if (w->stgs->verbose) {
      scs_printf("WARN: aa_init returned NULL, no acceleration applied.\n");
    }
  }
  return w;
}

static scs_int update_work(const ScsData *d, ScsWork *w,
                           const ScsSolution *sol) {
  /* before normalization */
  scs_int n = d->n;
  scs_int m = d->m;

  w->nm_b = SCS(norm)(d->b, m);
  w->nm_c = SCS(norm)(d->c, n);
  memcpy(w->b, d->b, d->m * sizeof(scs_float));
  memcpy(w->c, d->c, d->n * sizeof(scs_float));

#if EXTRA_VERBOSE > 0
  SCS(print_array)(w->b, m, "b");
  scs_printf("pre-normalized norm b = %4f\n", SCS(norm)(w->b, m));
  SCS(print_array)(w->c, n, "c");
  scs_printf("pre-normalized norm c = %4f\n", SCS(norm)(w->c, n));
#endif
  if (w->stgs->normalize) {
    SCS(normalize_b_c)(w);
#if EXTRA_VERBOSE > 0
    SCS(print_array)(w->b, m, "bn");
    scs_printf("sc_b = %4f\n", w->sc_b);
    scs_printf("post-normalized norm b = %4f\n", SCS(norm)(w->b, m));
    SCS(print_array)(w->c, n, "cn");
    scs_printf("sc_c = %4f\n", w->sc_c);
    scs_printf("post-normalized norm c = %4f\n", SCS(norm)(w->c, n));
#endif
  }
  if (w->stgs->warm_start) {
    warm_start_vars(w, sol);
  } else {
    cold_start_vars(w);
  }
  memcpy(w->h, w->c, n * sizeof(scs_float));
  memcpy(&(w->h[n]), w->b, m * sizeof(scs_float));
  memcpy(w->g, w->h, (n + m) * sizeof(scs_float));
  SCS(solve_lin_sys)(w->A, w->stgs, w->p, w->g, SCS_NULL, -1);
  SCS(scale_array)(&(w->g[n]), -1, m);
  w->g_th = SCS(dot)(w->h, w->g, n + m);
  return 0;
}

static scs_float iterate_norm_diff(ScsWork *w) {
  scs_int l = w->m + w->n + 1;
  scs_float u_norm_difference = SCS(norm_diff)(w->u, w->u_prev, l);
  scs_float v_norm_difference = SCS(norm_diff)(w->v, w->v_prev, l);
  scs_float norm = SQRTF(SCS(norm_sq)(w->u, l) + SCS(norm_sq)(w->v, l));
  scs_float norm_diff = SQRTF(u_norm_difference * u_norm_difference +
                              v_norm_difference * v_norm_difference);
  return norm_diff / norm;
}

static void update_best_iterate(ScsWork *w, ScsResiduals *r) {
  scs_float max_residual = get_max_residual(r);
  if (w->best_max_residual > max_residual) {
    w->best_max_residual = max_residual;
    memcpy(w->u_best, w->u, (w->m + w->n + 1) * sizeof(scs_float));
    memcpy(w->v_best, w->v, (w->m + w->n + 1) * sizeof(scs_float));
  }
}

scs_int SCS(solve)(ScsWork *w, const ScsData *d, const ScsCone *k,
                   ScsSolution *sol, ScsInfo *info) {
  scs_int i;
  SCS(timer) solve_timer, accel_timer;
  scs_float total_accel_time = 0.0, total_norm;
  ScsResiduals r;
  scs_int l = w->m + w->n + 1;
  if (!d || !k || !sol || !info || !w || !d->b || !d->c) {
    scs_printf("ERROR: SCS_NULL input\n");
    return SCS_FAILED;
  }
  /* initialize ctrl-c support */
  scs_start_interrupt_listener();
  SCS(tic)(&solve_timer);
  info->status_val = SCS_UNFINISHED; /* not yet converged */
  r.last_iter = -1;
  update_work(d, w, sol);

  if (w->stgs->verbose) {
    print_header(w, k);
  }
  /* scs: */
  for (i = 0; i < w->stgs->max_iters; ++i) {
    /* accelerate here so that last step always projection onto cone */
    /* this ensures the returned iterates always satisfy conic constraints */
    /* this relies on the fact that u and v are contiguous in memory */
    SCS(tic)(&accel_timer);
    if (i > 0 && aa_apply(w->u, w->u_prev, w->accel) != 0) {
      /*
      return failure(w, w->m, w->n, sol, info, SCS_FAILED,
          "error in accelerate", "Failure");
      */
    }
    total_accel_time += SCS(tocq)(&accel_timer);

    /* scs is homogeneous so scale the iterates to keep norm reasonable */
    total_norm = SQRTF(SCS(norm_sq)(w->u, l) + SCS(norm_sq)(w->v, l));
    SCS(scale_array)(w->u, SQRTF((scs_float)l) * ITERATE_NORM / total_norm, l);
    SCS(scale_array)(w->v, SQRTF((scs_float)l) * ITERATE_NORM / total_norm, l);

    memcpy(w->u_prev, w->u, l * sizeof(scs_float));
    memcpy(w->v_prev, w->v, l * sizeof(scs_float));

    if (project_lin_sys(w, i) < 0) {
      return failure(w, w->m, w->n, sol, info, SCS_FAILED,
                     "error in project_lin_sys", "Failure");
    }
    if (project_cones(w, k, i) < 0) {
      return failure(w, w->m, w->n, sol, info, SCS_FAILED,
                     "error in project_cones", "Failure");
    }

    update_dual_vars(w);

    if (scs_is_interrupted()) {
      return failure(w, w->m, w->n, sol, info, SCS_SIGINT, "Interrupted",
                     "Interrupted");
    }
    if (i % CONVERGED_INTERVAL == 0 || iterate_norm_diff(w) < 1e-10) {
      calc_residuals(w, &r, i);
      if ((info->status_val = has_converged(w, &r, i)) != 0) {
        break;
      }
      update_best_iterate(w, &r);
    }

    if (w->stgs->verbose && i % PRINT_INTERVAL == 0) {
      calc_residuals(w, &r, i);
      update_best_iterate(w, &r);
      print_summary(w, i, &r, &solve_timer);
    }
  }
  if (w->stgs->verbose) {
    calc_residuals(w, &r, i);
    print_summary(w, i, &r, &solve_timer);
  }
  /* populate solution vectors (unnormalized) and info */
  get_solution(w, sol, info, &r, i);
  info->solve_time = SCS(tocq)(&solve_timer);

  if (w->stgs->verbose) {
    print_footer(d, k, sol, w, info, total_accel_time);
  }
  scs_end_interrupt_listener();
  return info->status_val;
}

void SCS(finish)(ScsWork *w) {
  if (w) {
    SCS(finish_cone)(w->cone_work);
    if (w->stgs && w->stgs->normalize) {
#ifndef COPYAMATRIX
      SCS(un_normalize_a)(w->A, w->stgs, w->scal);
#else
      SCS(free_a_matrix)(w->A);
#endif
    }
    if (w->p) {
      SCS(free_lin_sys_work)(w->p);
    }
    if (w->accel) {
      aa_finish(w->accel);
    }
    free_work(w);
  }
}

ScsWork *SCS(init)(const ScsData *d, const ScsCone *k, ScsInfo *info) {
#if EXTRA_VERBOSE > 1
  SCS(tic)(&global_timer);
#endif
  ScsWork *w;
  SCS(timer) init_timer;
  scs_start_interrupt_listener();
  if (!d || !k || !info) {
    scs_printf("ERROR: Missing ScsData, ScsCone or ScsInfo input\n");
    return SCS_NULL;
  }
#if EXTRA_VERBOSE > 0
  SCS(print_data)(d);
  SCS(print_cone_data)(k);
#endif
#ifndef NOVALIDATE
  if (validate(d, k) < 0) {
    scs_printf("ERROR: Validation returned failure\n");
    return SCS_NULL;
  }
#endif
  SCS(tic)(&init_timer);
  if (d->stgs->write_data_filename) {
    SCS(write_data)(d, k);
  }
  w = init_work(d, k);
  info->setup_time = SCS(tocq)(&init_timer);
  if (d->stgs->verbose) {
    scs_printf("Setup time: %1.2es\n", info->setup_time / 1e3);
  }
  scs_end_interrupt_listener();
  return w;
}

/* this just calls SCS(init), SCS(solve), and SCS(finish) */
scs_int scs(const ScsData *d, const ScsCone *k, ScsSolution *sol,
            ScsInfo *info) {
  scs_int status;
  ScsWork *w = SCS(init)(d, k, info);
#if EXTRA_VERBOSE > 0
  scs_printf("size of scs_int = %lu, size of scs_float = %lu\n",
             sizeof(scs_int), sizeof(scs_float));
#endif
  if (w) {
    SCS(solve)(w, d, k, sol, info);
    status = info->status_val;
  } else {
    status = failure(SCS_NULL, d ? d->m : -1, d ? d->n : -1, sol, info,
                     SCS_FAILED, "could not initialize work", "Failure");
  }
  SCS(finish)(w);
  return status;
}
