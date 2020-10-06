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
    "   obj   ", " kap/tau ", " time (s)",
};
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 7;
static const scs_int LINE_LEN = 66;

static scs_int scs_isnan(scs_float x) { return (x == NAN || x != x); }

static void free_work(ScsWork *w) {
  if (w) {
    scs_free(w->u);
    scs_free(w->u_best);
    scs_free(w->u_t);
    scs_free(w->u_prev);
    scs_free(w->v);
    scs_free(w->v_best);
    scs_free(w->v_prev);
    scs_free(w->rsk);
    scs_free(w->h);
    scs_free(w->g);
    scs_free(w->b);
    scs_free(w->c);
    scs_free(w->ls_ws);
    scs_free(w->px);
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
  char *lin_sys_method = SCS(get_lin_sys_method)(d->A, d->P, d->stgs);
#ifdef USE_LAPACK
  scs_int acceleration_lookback = stgs->acceleration_lookback;
#else
  scs_int acceleration_lookback = 0;
#endif
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf(
      "\n\t       SCS v%s - Splitting Conic Solver\n\t(c) Brendan "
      "O'Donoghue, Stanford University, 2012\n",
      SCS(version)());
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  scs_printf("problem:  variables n: %i, constraints m: %i\n", (int)d->n,
             (int)d->m);
  scs_printf("%s", cone_str);
  scs_free(cone_str);
  scs_printf(
      "settings: eps: %.2e, alpha: %.2f, max_iters: %i,\n"
      "\t  normalize: %i, scale: %.2f, rho_x: %.2e,\n"
      "\t  acceleration_lookback: %i, warm_start: %i\n",
      stgs->eps, stgs->alpha, (int)stgs->max_iters, (int)stgs->normalize,
      stgs->scale, stgs->rho_x, (int)acceleration_lookback,
      (int)stgs->warm_start);
  if (lin_sys_method) {
    scs_printf("%s", lin_sys_method);
    scs_free(lin_sys_method);
  }

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
        sol->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
      }
      SCS(scale_array)(sol->x, NAN, n);
    }
    if (m > 0) {
      if (!sol->y) {
        sol->y = (scs_float *)scs_calloc(m, sizeof(scs_float));
      }
      SCS(scale_array)(sol->y, NAN, m);
      if (!sol->s) {
        sol->s = (scs_float *)scs_calloc(m, sizeof(scs_float));
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
  scs_float ax_s_sbctau, scale_i, pres = 0, *pr = w->pr;
  *nm_axs = 0;
  memset(pr, 0, w->m * sizeof(scs_float));
  SCS(accum_by_a)(w->A, w->p, x, pr);
  SCS(add_scaled_array)(pr, s, w->m, 1.0); /* pr = Ax + s */
  for (i = 0; i < w->m; ++i) {
    scale_i = 1.;
    if (w->stgs->normalize) {
      scale_i = 1. / w->scal->D[i] / w->scal->dual_scale;
    }
    *nm_axs += (pr[i] * scale_i) * (pr[i] * scale_i); /* ||Ax + s|| */
    ax_s_sbctau = (pr[i] - w->b[i] * tau) * scale_i;
    pres += ax_s_sbctau * ax_s_sbctau;
  }
  *nm_axs = SQRTF(*nm_axs);
  return SQRTF(pres); /* SCS(norm)(Ax + s - b * tau) */
}

/* we assume w->px contains unnormalized Px when this is called */
static scs_float calc_dual_resid(ScsWork *w, const scs_float *x,
                                 const scs_float *y, const scs_float tau,
                                 scs_float *nm_a_ty) {
  scs_int i;
  scs_float px_aty_ctau, scale_i, dres = 0, *dr = w->dr, *px = w->px;
  *nm_a_ty = 0;
  memset(w->dr, 0, w->n * sizeof(scs_float));
  SCS(accum_by_atrans)(w->A, w->p, y, dr); /* dr = A'y */
  for (i = 0; i < w->n; ++i) {
    scale_i = 1.;
    if (w->stgs->normalize) {
      scale_i = 1. / w->scal->E[i] / w->scal->primal_scale;
    }
    *nm_a_ty += (dr[i] * scale_i) * (dr[i] * scale_i); /* ||A' y|| */
    px_aty_ctau = (px[i] + dr[i] + w->c[i] * tau) * scale_i;
    dres += px_aty_ctau * px_aty_ctau;
  }
  *nm_a_ty = SQRTF(*nm_a_ty);
  return SQRTF(dres); /* SCS(norm)(Px + A'y + c * tau) */
}

/* calculates un-normalized quantities */
static void calc_residuals(ScsWork *w, ScsResiduals *r, scs_int iter) {
  scs_float *x = w->u, *y = &(w->u[w->n]);
  scs_float nmpr_tau, nmdr_tau, nm_axs_tau, nm_a_ty_tau, ct_x, bt_y;
  scs_float xt_p_x = 0.;
  scs_int n = w->n, m = w->m;

  /* checks if the residuals are unchanged by checking iteration */
  if (r->last_iter == iter) {
    return;
  }
  r->last_iter = iter;

  memset(w->px, 0, w->n * sizeof(scs_float));
  /* fills w->px with UNnormalized P * x */
  if (w->P) {
    /* px = P * x */
    SCS(accum_by_p)(w->P, w->p, x, w->px);
    /* xt_p_x = x' P x , unnormalized */
    xt_p_x = SCS(dot)(w->px, x, w->n);
  }

  r->tau = ABS(w->u[n + m]);
  r->kap = ABS(w->rsk[n + m]);

  nmpr_tau = calc_primal_resid(w, x, &(w->rsk[w->n]), r->tau, &nm_axs_tau);
  nmdr_tau = calc_dual_resid(w, x, y, r->tau, &nm_a_ty_tau);

  r->bt_y_by_tau = SCS(dot)(y, w->b, m);
  r->ct_x_by_tau = SCS(dot)(x, w->c, n);

  if (w->stgs->normalize) {
    r->kap /= (w->scal->primal_scale * w->scal->dual_scale);
    r->bt_y_by_tau /=
        (w->scal->primal_scale * w->scal->dual_scale);
    r->ct_x_by_tau /=
        (w->scal->primal_scale * w->scal->dual_scale);
    xt_p_x /= (w->scal->primal_scale * w->scal->dual_scale);
  }

  r->res_infeas = NAN;
  if (r->bt_y_by_tau < 0) {
    r->res_infeas = w->nm_b * nm_a_ty_tau / -r->bt_y_by_tau;
  }

  r->res_unbdd = NAN;
  r->xt_p_x_ctau = NAN;
  if (r->ct_x_by_tau < 0) {
    /* sqrt(x'Px) / (c'x) */
    r->xt_p_x_ctau = xt_p_x / r->ct_x_by_tau / r->ct_x_by_tau;
    /* |c||Ax + s| / (c'x) */
    r->res_unbdd = w->nm_c * nm_axs_tau / -r->ct_x_by_tau;
  }

  bt_y = SAFEDIV_POS(r->bt_y_by_tau, r->tau);
  ct_x = SAFEDIV_POS(r->ct_x_by_tau, r->tau);
  xt_p_x = SAFEDIV_POS(xt_p_x, r->tau * r->tau);

  r->xt_p_x = xt_p_x;
  r->res_pri = SAFEDIV_POS(nmpr_tau / (1 + w->nm_b), r->tau);
  r->res_dual = SAFEDIV_POS(nmdr_tau / (1 + w->nm_c), r->tau);
  r->rel_gap =
      ABS(xt_p_x + ct_x + bt_y) / (1 + ABS(xt_p_x) + ABS(ct_x) + ABS(bt_y));
}

static void cold_start_vars(ScsWork *w) {
  scs_int l = w->n + w->m + 1;
  memset(w->u, 0, l * sizeof(scs_float));
  memset(w->v, 0, l * sizeof(scs_float));
  w->u[l - 1] = SQRTF((scs_float)l);
  w->v[l - 1] = SQRTF((scs_float)l);
}

/* utility function that scales first n entries in inner prod by rho_x */
/* and last m entries by 1 / scale, assumes length of array is n + m */
static scs_float dot_with_diag_scaling(ScsWork *w, const scs_float *x,
                                const scs_float *y) {
  scs_int i, n = w->n, len = w->n + w->m;
  scs_float ip = 0.0;
  for (i = 0; i < n; ++i) {
    ip += w->stgs->rho_x * x[i] * y[i];
  }
  for (i = n; i < len; ++i) {
    ip += x[i] * y[i] / w->stgs->scale;
  }
  return ip;
}

static scs_float root_plus(ScsWork *w, scs_float *p, scs_float *mu,
                           scs_float eta) {
  scs_float b, c, tau, a = w->root_plus_a;
  b = dot_with_diag_scaling(w, mu, w->g) - 2 * dot_with_diag_scaling(w, p, w->g) - eta;
  c = dot_with_diag_scaling(w, p, p) - dot_with_diag_scaling(w, p, mu);
  tau = (-b + SQRTF(MAX(b * b - 4 * a * c, 0.))) / (2 * a);
#if EXTRA_VERBOSE > 3
  scs_printf("root_plus: a: %g, b: %g, c: %g, eta: %g, tau: %g, tau no p: %g\n",
             a, b, c, eta, tau,
             MAX(0., (eta + SCS(dot)(p, w->h, w->m + w->n)) /
                         (1 + SCS(dot)(w->h, w->g, w->m + w->n))));
#endif
  return tau;
}

/*
 * Forms a warm start for the linear solve step using the output of the cone
 * projection step.
 */
static void compute_warm_start_lin_sys(ScsWork *w) {
  scs_int len = w->m + w->n;
  memcpy(w->ls_ws, w->u, len * sizeof(scs_float));
  SCS(add_scaled_array)(w->ls_ws, w->g, len, w->u[len]);
}

/* status < 0 indicates failure */
static scs_int project_lin_sys(ScsWork *w, scs_int iter) {
  scs_int n = w->n, m = w->m, l = n + m + 1, status;
  memcpy(w->u_t, w->v, l * sizeof(scs_float));
  SCS(scale_array)(w->u_t, w->stgs->rho_x, n);
  SCS(scale_array)(&(w->u_t[n]), -1. / w->stgs->scale, m);
  /* compute warm start, if used */
  compute_warm_start_lin_sys(w);
  /* now w->ls_ws contains warm start */
  status =
      SCS(solve_lin_sys)(w->A, w->P, w->stgs, w->p, w->u_t, w->ls_ws, iter);
  w->u_t[l - 1] = root_plus(w, w->u_t, w->v, w->v[l - 1]);
  SCS(add_scaled_array)(w->u_t, w->g, l - 1, -w->u_t[l - 1]);
  return status;
}

/* compute the [r;s;kappa] iterate */
/* rsk^{k+1} = u^{k+1} + v^k - 2 u_t^{k+1} */
/* uses Moreau decomposition to get projection onto dual cone */
/* since it depends on v^k MUST be called before update_dual_vars is done */
/* effect of w->stgs->alpha is cancelled out */
static void compute_rsk(ScsWork *w) {
  scs_int i, l = w->m + w->n + 1;
  /* r, should = 0 */
  for (i = 0; i < w->n; ++i) {
    w->rsk[i] = w->stgs->rho_x * (w->v[i] + w->u[i] - 2 * w->u_t[i]);
  }
#if EXTRA_VERBOSE > 4
  scs_printf("norm(r) = %1.3f\n", SCS(norm)(w->rsk, w->n));
#endif
  /* s */
  for (i = w->n; i < l-1; ++i) {
    w->rsk[i] = (w->v[i] + w->u[i] - 2 * w->u_t[i]) / w->stgs->scale;
  }
  /* kappa */
  w->rsk[l-1] = w->v[l-1] + w->u[l-1] - 2 * w->u_t[l-1];
}

static void update_dual_vars(ScsWork *w) {
  scs_int i, l = w->n + w->m + 1;
  scs_float a = w->stgs->alpha;
  /* compute and store [r;s;kappa] */
  compute_rsk(w);
  for (i = 0; i < l; ++i) {
    w->v[i] += a * (w->u[i] - w->u_t[i]);
  }
}

/* status < 0 indicates failure */
static scs_int project_cones(ScsWork *w, const ScsCone *k, scs_int iter) {
  scs_int i, n = w->n, l = w->n + w->m + 1, st;
  for (i = 0; i < l; ++i) {
    w->u[i] = 2 * w->u_t[i] - w->v[i];
  }
  /* u = [x;y;tau] */
  st = SCS(proj_dual_cone)(&(w->u[n]), k, w->cone_work, &(w->u_prev[n]), iter,
                           w->scal->D);
  if (w->u[l - 1] < 0.0) {
    w->u[l - 1] = 0.0;
  }
  return st;
}

static scs_int indeterminate(ScsWork *w, ScsSolution *sol, ScsInfo *info) {
  strcpy(info->status, "indeterminate");
  SCS(scale_array)(sol->x, NAN, w->n);
  SCS(scale_array)(sol->y, NAN, w->m);
  SCS(scale_array)(sol->s, NAN, w->m);
  return SCS_INDETERMINATE;
}

static void sety(ScsWork *w, ScsSolution *sol) {
  if (!sol->y) {
    sol->y = (scs_float *)scs_calloc(w->m, sizeof(scs_float));
  }
  memcpy(sol->y, &(w->u[w->n]), w->m * sizeof(scs_float));
}

/* s is contained in rsk */
static void sets(ScsWork *w, ScsSolution *sol) {
  if (!sol->s) {
    sol->s = (scs_float *)scs_calloc(w->m, sizeof(scs_float));
  }
  memcpy(sol->s, &(w->rsk[w->n]), w->m * sizeof(scs_float));
}

static void setx(ScsWork *w, ScsSolution *sol) {
  if (!sol->x) {
    sol->x = (scs_float *)scs_calloc(w->n, sizeof(scs_float));
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
    strcpy(info->status, "solved / inaccurate");
    return SCS_SOLVED_INACCURATE;
  }
  strcpy(info->status, "solved");
  return SCS_SOLVED;
}

static scs_int infeasible(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                          scs_float bt_y) {
  SCS(scale_array)(sol->y, -1 / bt_y, w->m);
  SCS(scale_array)(sol->x, NAN, w->n);
  SCS(scale_array)(sol->s, NAN, w->m);
  if (info->status_val == 0) {
    strcpy(info->status, "infeasible / inaccurate");
    return SCS_INFEASIBLE_INACCURATE;
  }
  strcpy(info->status, "infeasible");
  return SCS_INFEASIBLE;
}

static scs_int unbounded(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                         scs_float ct_x) {
  SCS(scale_array)(sol->x, -1 / ct_x, w->n);
  SCS(scale_array)(sol->s, -1 / ct_x, w->m);
  SCS(scale_array)(sol->y, NAN, w->m);
  if (info->status_val == 0) {
    strcpy(info->status, "unbounded / inaccurate");
    return SCS_UNBOUNDED_INACCURATE;
  }
  strcpy(info->status, "unbounded");
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
    info->xt_p_x = r->xt_p_x;
    info->pobj = r->xt_p_x / 2. + r->ct_x_by_tau / r->tau;
    info->dobj = -r->xt_p_x / 2. - r->bt_y_by_tau / r->tau;
  } else if (is_unbounded_status(info->status_val)) {
    info->rel_gap = NAN;
    info->res_pri = NAN;
    info->res_dual = NAN;
    info->xt_p_x = r->xt_p_x_ctau;
    info->pobj = -INFINITY;
    info->dobj = -INFINITY;
  } else if (is_infeasible_status(info->status_val)) {
    info->rel_gap = NAN;
    info->res_pri = NAN;
    info->res_dual = NAN;
    info->xt_p_x = NAN;
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
  scs_float pobj = r->xt_p_x / 2. + SAFEDIV_POS(r->ct_x_by_tau, r->tau);
  scs_printf("%*i|", (int)strlen(HEADER[0]), (int)i);
  scs_printf("%*.2e ", (int)HSPACE, r->res_pri);
  scs_printf("%*.2e ", (int)HSPACE, r->res_dual);
  scs_printf("%*.2e ", (int)HSPACE, r->rel_gap);
  scs_printf("%*.2e ", (int)HSPACE, pobj);
  scs_printf("%*.2e ", (int)HSPACE, SAFEDIV_POS(r->kap, r->tau));
  scs_printf("%*.2e ", (int)HSPACE, SCS(tocq)(solve_timer) / 1e3);
  scs_printf("\n");

#if EXTRA_VERBOSE > 0
  scs_printf("Norm u = %4f, ", SCS(norm)(w->u, w->n + w->m + 1));
  scs_printf("Norm u_t = %4f, ", SCS(norm)(w->u_t, w->n + w->m + 1));
  scs_printf("Norm v = %4f, ", SCS(norm)(w->v, w->n + w->m + 1));
  scs_printf("tau = %4f, ", w->u[w->n + w->m]);
  scs_printf("kappa = %4f, ", w->rsk[w->n + w->m]);
  scs_printf("|u - u_prev| = %1.2e, ",
             SCS(norm_diff)(w->u, w->u_prev, w->n + w->m + 1));
  scs_printf("|u - u_t| = %1.2e, ",
             SCS(norm_diff)(w->u, w->u_t, w->n + w->m + 1));
  scs_printf("res_infeas = %1.2e, ", r->res_infeas);
  scs_printf("res_unbdd = %1.2e, ", r->res_unbdd);
  scs_printf("xt_p_x_ctau = %1.2e\n", r->xt_p_x_ctau);
#endif

#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

static void print_header(ScsWork *w, const ScsCone *k) {
  scs_int i;
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

// XXX norm inf change?
static scs_float get_dual_cone_dist(const scs_float *y, const ScsCone *k,
                                    ScsConeWork *c, scs_int m) {
  scs_float dist;
  scs_float *t = (scs_float *)scs_calloc(m, sizeof(scs_float));
  memcpy(t, y, m * sizeof(scs_float));
  SCS(proj_dual_cone)(t, k, c, SCS_NULL, -1, SCS_NULL);
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
// XXX norm inf change?
static scs_float get_pri_cone_dist(const scs_float *s, const ScsCone *k,
                                   ScsConeWork *c, scs_int m) {
  scs_float dist;
  scs_float *t = (scs_float *)scs_calloc(m, sizeof(scs_float));
  memcpy(t, s, m * sizeof(scs_float));
  SCS(scale_array)(t, -1.0, m);
  SCS(proj_dual_cone)(t, k, c, SCS_NULL, -1, SCS_NULL);
  dist = SCS(norm_inf)(t, m); /* ||s - Pi_c(s)|| = ||Pi_c*(-s)|| */
#if EXTRA_VERBOSE > 0
  SCS(print_array)(s, m, "s");
  SCS(print_array)(t, m, "(s - proj_s)");
  scs_printf("dist = %4f\n", dist);
#endif
  scs_free(t);
  return dist;
}

static void print_footer(const ScsData *d, const ScsCone *k, ScsSolution *sol,
                         ScsWork *w, ScsInfo *info,
                         scs_float total_lin_sys_time,
                         scs_float total_cone_time,
                         scs_float total_accel_time) {
  scs_int i;
  char *lin_sys_str = SCS(get_lin_sys_summary)(w->p, info);

  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  if (info->iter == w->stgs->max_iters) {
    scs_printf("hit max_iters, returning best iterate\n");
  }
  scs_printf("status:  %s\n", info->status);
  scs_printf("timings: total: %1.2es = setup: %1.2es + solve: %1.2es\n",
             (info->setup_time + info->solve_time) / 1e3,
             info->setup_time / 1e3, info->solve_time / 1e3);
  scs_printf("\t lin-sys: %1.2es, cones: %1.2es, accel: %1.2es\n",
             total_lin_sys_time / 1e3, total_cone_time / 1e3,
             total_accel_time / 1e3);
  scs_printf("%s", lin_sys_str);
  scs_free(lin_sys_str);

  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");

  if (is_infeasible_status(info->status_val)) {
    scs_printf("cone: dist(y, K*) = %.4e\n",
               get_dual_cone_dist(sol->y, k, w->cone_work, d->m));
    scs_printf("cert: |A'y|_2*|b|_2 = %.4e\n", info->res_infeas);
    scs_printf("      b'y = %.4f\n", SCS(dot)(d->b, sol->y, d->m));
  } else if (is_unbounded_status(info->status_val)) {
    scs_printf("cone: dist(s, K) = %.4e\n",
               get_pri_cone_dist(sol->s, k, w->cone_work, d->m));
    scs_printf("cert: |Ax+s|_2*|c|_2 = %.4e\n", info->res_unbdd);
    scs_printf("      (x'Px)^(1/2) = %.4e\n", SQRTF(MAX(info->xt_p_x, 0.)));
    scs_printf("      c'x = %.4f\n", SCS(dot)(d->c, sol->x, d->n));
  } else {
    scs_printf("cones: dist(s, K) = %.4e, dist(y, K*) = %.4e\n",
               get_pri_cone_dist(sol->s, k, w->cone_work, d->m),
               get_dual_cone_dist(sol->y, k, w->cone_work, d->m));
    scs_printf("comp slack: s'y/|s||y| = %.4e\n",
               SCS(dot)(sol->s, sol->y, d->m) /
               MAX(1e-9, SCS(norm)(sol->s, d->m)) /
               MAX(1e-9, SCS(norm)(sol->y, d->m)));
    scs_printf("primal res: |Ax+s-b|_2/(1+|b|_2) = %.4e\n", info->res_pri);
    scs_printf("dual res: |Px+A'y+c|_2/(1+|c|_2) = %.4e\n", info->res_dual);
    scs_printf("rel gap: |x'Px+c'x+b'y|/(1+|x'Px|+|c'x|+|b'y|) = %.4e\n",
               info->rel_gap);
  }
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  scs_printf("optimal objective = %.6f\n", info->pobj);
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
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
  if (isless(r->res_unbdd, eps) &&
      isless(SQRTF(MAX(r->xt_p_x_ctau, 0.)), eps) && iter > 0) {
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
    /* scs_printf("WARN: m less than n, problem likely degenerate\n"); */
    /* return -1; */
  }
  if (SCS(validate_lin_sys)(d->A, d->P) < 0) {
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
  w->last_scale_update_iter = 0;
  w->log_scale_factor_mean = 0.;
  w->best_max_residual = INFINITY;
  /* allocate workspace: */
  w->u = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->u_best = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->u_t = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->u_prev = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v_prev = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v_best = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->rsk = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->h = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->g = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->ls_ws = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->px = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  w->pr = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  w->dr = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  w->b = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  w->c = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  if (!w->u || !w->u_t || !w->u_prev || !w->h || !w->g || !w->pr || !w->dr ||
      !w->b || !w->c) {
    scs_printf("ERROR: work memory allocation failure\n");
    return SCS_NULL;
  }
  if (!(w->cone_work = SCS(init_cone)(k))) {
    scs_printf("ERROR: init_cone failure\n");
    return SCS_NULL;
  }
  w->A = d->A;
  w->P = d->P;
  if (w->stgs->normalize) {
#ifdef COPYAMATRIX
    if (!SCS(copy_matrix)(&(w->A), d->A)) {
      scs_printf("ERROR: copy A matrix failed\n");
      return SCS_NULL;
    }
    if (w->P && !SCS(copy_matrix)(&(w->P), d->P)) {
      scs_printf("ERROR: copy P matrix failed\n");
      return SCS_NULL;
    }
#endif
    w->scal = (ScsScaling *)scs_calloc(1, sizeof(ScsScaling));
    SCS(normalize)(w->A, w->P, k, w->scal, w->cone_work);
#if EXTRA_VERBOSE > 0
    SCS(print_array)(w->scal->D, d->m, "D");
    scs_printf("norm(D) = %4f\n", SCS(norm)(w->scal->D, d->m));
    SCS(print_array)(w->scal->E, d->n, "E");
    scs_printf("norm(E) = %4f\n", SCS(norm)(w->scal->E, d->n));
#endif
  } else {
    w->scal = SCS_NULL;
  }
  if (!(w->p = SCS(init_lin_sys_work)(w->A, w->P, w->stgs))) {
    scs_printf("ERROR: init_lin_sys_work failure\n");
    return SCS_NULL;
  }
  /* hack: negative acceleration_lookback interpreted as type-I */
  if (!(w->accel = aa_init(l, ABS(w->stgs->acceleration_lookback),
                           w->stgs->acceleration_lookback < 0, ETA))) {
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
    scs_printf("dual scale= %4f\n", w->scal->dual_scale);
    scs_printf("post-normalized norm b = %4f\n", SCS(norm)(w->b, m));
    SCS(print_array)(w->c, n, "cn");
    scs_printf("primal scale= %4f\n", w->scal->primal_scale);
    scs_printf("post-normalized norm c = %4f\n", SCS(norm)(w->c, n));
#endif
  }
  if (w->stgs->warm_start) {
    warm_start_vars(w, sol);
  } else {
    cold_start_vars(w);
  }

  /* h = [c;b] */
  memcpy(w->h, w->c, n * sizeof(scs_float));
  memcpy(&(w->h[n]), w->b, m * sizeof(scs_float));

  /* g = (I + M)^{-1} h */
  memcpy(w->g, w->h, (n + m) * sizeof(scs_float));
  SCS(scale_array)(&(w->g[n]), -1., m);
  SCS(solve_lin_sys)(w->A, w->P, w->stgs, w->p, w->g, SCS_NULL, -1);
  w->root_plus_a = 1 + dot_with_diag_scaling(w, w->g, w->g);
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

static void maybe_update_scale(ScsWork *w, ScsResiduals *r, scs_int iter) {
  scs_float factor;
  // TODO XXX
  scs_int iters_since_last_update = (iter - w->last_scale_update_iter) / 20;
  w->log_scale_factor_mean *= iters_since_last_update;
  /* higher scale makes res_pri go down faster, so increase is res_pri larger */
  w->log_scale_factor_mean += log(r->res_pri / r->res_dual);
  w->log_scale_factor_mean /= (iters_since_last_update + 1);
  //scs_printf("ratio %4f\n", log(r->res_pri / r->res_dual));
  //scs_printf("log_scale_factor_mean %4f\n", w->log_scale_factor_mean);
  //scs_printf("iters_since_last_update %i\n", iters_since_last_update);
  /* XXX: bound this ? */
  factor = SQRTF(exp(w->log_scale_factor_mean));
  //scs_printf("factor %4f\n", factor);
  if (SCS(should_update_scale(factor, iter))) {
    w->stgs->scale *= factor;
    scs_printf("new scale %4f\n", w->stgs->scale);
    w->log_scale_factor_mean = 0;
    w->last_scale_update_iter = iter;

    SCS(update_linsys_scale)(w->A, w->P, w->stgs, w->p);
    /* g = (I + M)^{-1} h */
    memcpy(w->g, w->h, (w->n + w->m) * sizeof(scs_float));
    SCS(scale_array)(&(w->g[w->n]), -1., w->m);
    SCS(solve_lin_sys)(w->A, w->P, w->stgs, w->p, w->g, SCS_NULL, -1);
    w->root_plus_a = 1 + dot_with_diag_scaling(w, w->g, w->g);

    /* XXX reset aa? */
    /* XXX update v somehow? */
    return;
  }
}

scs_int SCS(solve)(ScsWork *w, const ScsData *d, const ScsCone *k,
                   ScsSolution *sol, ScsInfo *info) {
  scs_int i;
  scs_float v_norm;
  SCS(timer) solve_timer, lin_sys_timer, cone_timer, accel_timer;
  scs_float total_accel_time = 0.0, total_cone_time = 0.0,
            total_lin_sys_time = 0.0;
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

  /* SCS */
  for (i = 0; i < w->stgs->max_iters; ++i) {
    /* scs is homogeneous so scale the iterate to keep norm reasonable */
    v_norm = SCS(norm)(w->v, l);
    SCS(scale_array)(w->v, SQRTF((scs_float)l) * ITERATE_NORM / v_norm, l);

    /* XXX rm this? */
    memcpy(w->u_prev, w->u, l * sizeof(scs_float));
    memcpy(w->v_prev, w->v, l * sizeof(scs_float));

    SCS(tic)(&lin_sys_timer);
    if (project_lin_sys(w, i) < 0) {
      return failure(w, w->m, w->n, sol, info, SCS_FAILED,
                     "error in project_lin_sys", "Failure");
    }
    total_lin_sys_time += SCS(tocq)(&lin_sys_timer);

    SCS(tic)(&cone_timer);
    if (project_cones(w, k, i) < 0) {
      return failure(w, w->m, w->n, sol, info, SCS_FAILED,
                     "error in project_cones", "Failure");
    }
    total_cone_time += SCS(tocq)(&cone_timer);

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

    /* Finally apply any acceleration */
    SCS(tic)(&accel_timer);
    if (aa_apply(w->v, w->v_prev, w->accel) != 0) {
      /*
      return failure(w, w->m, w->n, sol, info, SCS_FAILED,
          "error in accelerate", "Failure");
      */
    }
    total_accel_time += SCS(tocq)(&accel_timer);

    /* if residuals are fresh then maybe compute new scale */
    if (i == r.last_iter) {
      maybe_update_scale(w, &r, i);
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
    print_footer(d, k, sol, w, info, total_lin_sys_time, total_cone_time,
                 total_accel_time);
  }

  scs_end_interrupt_listener();
  return info->status_val;
}

void SCS(finish)(ScsWork *w) {
  if (w) {
    SCS(finish_cone)(w->cone_work);
    if (w->stgs && w->stgs->normalize) {
#ifndef COPYAMATRIX
      SCS(un_normalize)(w->A, w->P, w->scal);
#else
      SCS(free_scs_matrix)(w->A);
      SCS(free_scs_matrix)(w->P);
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
