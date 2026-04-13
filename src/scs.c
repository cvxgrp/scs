/*
 * Main SCS solver implementation.
 *
 * This file contains the ADMM iteration loop, workspace management,
 * residual computation, convergence checking, and the public API
 * functions (scs_init, scs_solve, scs_update, scs_finish, scs).
 */

#include "scs.h"
#include "aa.h"
#include "cones.h"
#include "ctrlc.h"
#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "normalize.h"
#include "rw.h"
#include "scs_matrix.h"
#include "scs_work.h"
#include "util.h"

#include <string.h>

#ifdef SCS_MKL
#define MKL_INTERFACE_LP64 0
#define MKL_INTERFACE_ILP64 1
int MKL_Set_Interface_Layer(int);

static const char *scs_mkl_interface_name(int layer) {
  switch (layer) {
  case MKL_INTERFACE_LP64:
    return "LP64";
  case MKL_INTERFACE_ILP64:
    return "ILP64";
  default:
    return "unknown";
  }
}

static scs_int scs_init_mkl_runtime(void) {
  /* Enforce the correct MKL interface layer for BLAS/LAPACK calls before any
   * common SCS code reaches MKL-backed BLAS routines.
   *
   * The interface layer must match what we linked against:
   *   BLAS64 defined   -> ILP64 (64-bit BLAS integers, mkl-dynamic-ilp64-*)
   *   BLAS64 undefined -> LP64  (32-bit BLAS integers, mkl-dynamic-lp64-*)
   *
   * If another library in the process set the wrong layer, our BLAS calls
   * would silently receive the wrong integer width, causing memory corruption.
   *
   * This only protects the BLAS layer. The pardiso_64 entry point is
   * unaffected by the interface layer — it always uses 64-bit integers. */
  {
#ifdef BLAS64
    int expected = MKL_INTERFACE_ILP64;
#else
    int expected = MKL_INTERFACE_LP64;
#endif
    int actual = MKL_Set_Interface_Layer(expected);
    if (actual != expected) {
      scs_printf("MKL interface layer mismatch: expected %s, but MKL is using "
                 "%s (%d). Another library in this process likely "
                 "initialized MKL with an incompatible LP64/ILP64 setting.\n",
                 scs_mkl_interface_name(expected),
                 scs_mkl_interface_name(actual), actual);
      return -1;
    }
  }
  return 0;
}
#endif

/* ======================= Forward Declarations ====================== */

static void print_init_header(const ScsData *d, const ScsCone *k,
                              const ScsSettings *stgs);
static void print_header(ScsWork *w, const ScsCone *k);
static void print_summary(ScsWork *w, scs_int i, SCS(timer) * solve_timer);
static void print_footer(ScsInfo *info);
static void free_residuals(ScsResiduals *r);
static ScsResiduals *init_residuals(const ScsData *d);
static void populate_on_failure(scs_int m, scs_int n, ScsSolution *sol,
                                ScsInfo *info, scs_int status_val,
                                const char *msg);
static scs_int failure(ScsWork *w, scs_int m, scs_int n, ScsSolution *sol,
                       ScsInfo *info, scs_int stint, const char *msg,
                       const char *ststr);
static void compute_residuals(ScsResiduals *r, scs_int m, scs_int n,
                              scs_float pd);
static void unnormalize_residuals(ScsWork *w);
static void populate_residual_struct(ScsWork *w, scs_int iter);
static scs_int has_converged(ScsWork *w, scs_int iter);
static void warm_start_vars(ScsWork *w, ScsSolution *sol);
static void cold_start_vars(ScsWork *w);
static scs_float root_plus(ScsWork *w, scs_float *p, scs_float *mu,
                           scs_float eta);
static scs_int project_lin_sys(ScsWork *w, scs_int iter);
static void compute_rsk(ScsWork *w);
static void update_dual_vars(ScsWork *w);
static scs_int project_cones(ScsWork *w, const ScsCone *k, scs_int iter);
static void finalize(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                     scs_int iter);
static void set_diag_r(ScsWork *w);
static ScsWork *init_work(const ScsData *d, const ScsCone *k,
                          const ScsSettings *stgs);
static void update_work_cache(ScsWork *w);
static void reset_tracking(ScsWork *w);
static scs_int update_work(ScsWork *w, ScsSolution *sol);
static scs_int update_scale(ScsWork *w, const ScsCone *k, scs_int iter);
static inline void normalize_v(scs_float *v, scs_int len);

/* ======================== Printing / Output ======================== */

static const char *HEADER[] = {
    " iter ",    " pri res ", " dua res ", "   gap   ",
    "   obj   ", "  scale  ", " time (s)",
};
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 7;
static const scs_int LINE_LEN = 66;

static void print_init_header(const ScsData *d, const ScsCone *k,
                              const ScsSettings *stgs) {
  scs_int i;
  char *cone_str = SCS(get_cone_header)(k);
  const char *lin_sys_method = scs_get_lin_sys_method();
#ifdef USE_LAPACK
  scs_int acceleration_lookback = stgs->acceleration_lookback;
  scs_int acceleration_interval = stgs->acceleration_interval;
#else
  scs_int acceleration_lookback = 0;
  scs_int acceleration_interval = 0;
#endif
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n\t       SCS v%s - Splitting Conic Solver\n\t(c) Brendan "
             "O'Donoghue, Stanford University, 2012\n",
             scs_version());
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  scs_printf("problem:  variables n: %i, constraints m: %i\n", (int)d->n,
             (int)d->m);
  if (cone_str) {
    scs_printf("%s", cone_str);
    scs_free(cone_str);
  } else {
    scs_printf("cones: <unavailable>\n");
  }
  scs_printf("settings: eps_abs: %.1e, eps_rel: %.1e, eps_infeas: %.1e\n"
             "\t  alpha: %.2f, scale: %.2e, adaptive_scale: %i\n"
             "\t  max_iters: %i, normalize: %i, rho_x: %.2e\n",
             stgs->eps_abs, stgs->eps_rel, stgs->eps_infeas, stgs->alpha,
             stgs->scale, (int)stgs->adaptive_scale, (int)stgs->max_iters,
             (int)stgs->normalize, stgs->rho_x);
  if (stgs->acceleration_lookback != 0) {
    scs_printf("\t  acceleration_lookback: %i, acceleration_interval: %i\n",
               (int)acceleration_lookback, (int)acceleration_interval);
  }
  if (stgs->time_limit_secs) {
    scs_printf("\t  time_limit_secs: %.2e\n", stgs->time_limit_secs);
  }
#ifdef _OPENMP
  scs_printf("\t  compiled with openmp parallelization enabled\n");
#endif
  if (lin_sys_method) {
    scs_printf("lin-sys:  %s\n\t  nnz(A): %li, nnz(P): %li\n", lin_sys_method,
               (long)d->A->p[d->A->n], d->P ? (long)d->P->p[d->P->n] : 0l);
  }

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

static void print_summary(ScsWork *w, scs_int i, SCS(timer) * solve_timer) {
  ScsResiduals *r = w->r_orig;
  scs_printf("%*i|", (int)strlen(HEADER[0]), (int)i);
  scs_printf("%*.2e ", (int)HSPACE, r->res_pri);
  scs_printf("%*.2e ", (int)HSPACE, r->res_dual);
  scs_printf("%*.2e ", (int)HSPACE, r->gap);
  /* report mid point of primal and dual objective values */
  scs_printf("%*.2e ", (int)HSPACE, 0.5 * (r->pobj + r->dobj));
  scs_printf("%*.2e ", (int)HSPACE, w->stgs->scale);
  /* Report TOTAL time, including setup */
  scs_printf("%*.2e ", (int)HSPACE,
             (SCS(tocq)(solve_timer) + w->setup_time) / 1e3);
  scs_printf("\n");

#if VERBOSITY > 0
  scs_printf("Norm u = %1.6e, ", SCS(norm_2)(w->u, w->d->n + w->d->m + 1));
  scs_printf("Norm u_t = %1.6e, ", SCS(norm_2)(w->u_t, w->d->n + w->d->m + 1));
  scs_printf("Norm v = %1.6e, ", SCS(norm_2)(w->v, w->d->n + w->d->m + 1));
  scs_printf("Norm rsk = %1.6e, ", SCS(norm_2)(w->rsk, w->d->n + w->d->m + 1));
  scs_printf("Norm x = %1.6e, ", SCS(norm_2)(w->xys_orig->x, w->d->n));
  scs_printf("Norm y = %1.6e, ", SCS(norm_2)(w->xys_orig->y, w->d->m));
  scs_printf("Norm s = %1.6e, ", SCS(norm_2)(w->xys_orig->s, w->d->m));
  scs_printf("Norm |Ax + s| = %1.6e, ", SCS(norm_2)(r->ax_s, w->d->m));
  scs_printf("tau = %1.6e, ", w->u[w->d->n + w->d->m]);
  scs_printf("kappa = %1.6e, ", w->rsk[w->d->n + w->d->m]);
  scs_printf("|u - u_t| = %1.6e, ",
             SCS(norm_diff)(w->u, w->u_t, w->d->n + w->d->m + 1));
  scs_printf("res_infeas = %1.6e, ", r->res_infeas);
  scs_printf("res_unbdd_a = %1.6e, ", r->res_unbdd_a);
  scs_printf("res_unbdd_p = %1.6e, ", r->res_unbdd_p);
  scs_printf("ctx_tau = %1.6e, ", r->ctx_tau);
  scs_printf("bty_tau = %1.2e\n", r->bty_tau);
#endif

#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

static void print_footer(ScsInfo *info) {
  scs_int i;

  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  scs_printf("status:  %s\n", info->status);
  scs_printf("timings: total: %1.2es = setup: %1.2es + solve: %1.2es\n",
             (info->setup_time + info->solve_time) / 1e3,
             info->setup_time / 1e3, info->solve_time / 1e3);
  scs_printf("\t lin-sys: %1.2es, cones: %1.2es, accel: %1.2es\n",
             info->lin_sys_time / 1e3, info->cone_time / 1e3,
             info->accel_time / 1e3);

  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
  /* report mid point of primal and dual objective values */
  scs_printf("objective = %.6f", 0.5 * (info->pobj + info->dobj));
  switch (info->status_val) {
  case SCS_SOLVED_INACCURATE:
  case SCS_UNBOUNDED_INACCURATE:
  case SCS_INFEASIBLE_INACCURATE:
    scs_printf(" (inaccurate)");
    /* fallthrough */
  default:
    scs_printf("\n");
  }
  for (i = 0; i < LINE_LEN; ++i) {
    scs_printf("-");
  }
  scs_printf("\n");
#ifdef MATLAB_MEX_FILE
  mexEvalString("drawnow;");
#endif
}

/* ======================== Memory Management ======================== */

static void free_residuals(ScsResiduals *r) {
  if (r) {
    scs_free(r->ax);
    scs_free(r->ax_s);
    scs_free(r->px);
    scs_free(r->aty);
    scs_free(r->ax_s_btau);
    scs_free(r->px_aty_ctau);
    scs_free(r);
  }
}

static ScsResiduals *init_residuals(const ScsData *d) {
  ScsResiduals *r = (ScsResiduals *)scs_calloc(1, sizeof(ScsResiduals));
  if (!r)
    return SCS_NULL;
  r->ax = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  r->ax_s = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  r->ax_s_btau = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  r->px = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  r->aty = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  r->px_aty_ctau = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  return r;
}

/* ==================== Error / Failure Handling ===================== */

static void populate_on_failure(scs_int m, scs_int n, ScsSolution *sol,
                                ScsInfo *info, scs_int status_val,
                                const char *msg) {
  if (info) {
    info->gap = NAN;
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

/* ========================= Validation ============================== */

#if NO_VALIDATE == 0
static scs_int validate(const ScsData *d, const ScsCone *k,
                        const ScsSettings *stgs) {
  if (d->m <= 0 || d->n <= 0) {
    scs_printf("m and n must both be greater than 0; m = %li, n = %li\n",
               (long)d->m, (long)d->n);
    return -1;
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
  if (stgs->eps_abs < 0) {
    scs_printf("eps_abs tolerance must be positive\n");
    return -1;
  }
  if (stgs->eps_rel < 0) {
    scs_printf("eps_rel tolerance must be positive\n");
    return -1;
  }
  if (stgs->eps_infeas < 0) {
    scs_printf("eps_infeas tolerance must be positive\n");
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
  if (stgs->acceleration_interval <= 0) {
    scs_printf("acceleration_interval must be positive (10 works well).\n");
    return -1;
  }
  return 0;
}
#endif

/* ==================== Residual Computation ========================= */

/* pd = primal_scale * dual_scale = sigma^2 from the Ruiz equilibration.
 * Pass pd = 1 when operating on normalised residuals (no amplification).
 *
 * After unnormalization bty_tau and ctx_tau are divided by pd, so
 * floating-point noise is amplified by 1/pd. INFEAS_NEGATIVITY_TOL is
 * calibrated for normalised scale, so the correct threshold in the
 * un-normalised space is INFEAS_NEGATIVITY_TOL/pd (issue #350). */
static void compute_residuals(ScsResiduals *r, scs_int m, scs_int n,
                              scs_float pd) {
  scs_float nm_ax_s, nm_px, nm_aty;
  scs_float nm_ax_s_btau = NORM(r->ax_s_btau, m);
  scs_float nm_px_aty_ctau = NORM(r->px_aty_ctau, n);
  scs_float tol = INFEAS_NEGATIVITY_TOL / pd;

  r->res_pri = SAFEDIV_POS(nm_ax_s_btau, r->tau);
  r->res_dual = SAFEDIV_POS(nm_px_aty_ctau, r->tau);
  r->res_unbdd_a = NAN;
  r->res_unbdd_p = NAN;
  r->res_infeas = NAN;
  if (r->ctx_tau < -tol) {
    nm_ax_s = NORM(r->ax_s, m);
    nm_px = NORM(r->px, n);
    r->res_unbdd_a = SAFEDIV_POS(nm_ax_s, -r->ctx_tau);
    r->res_unbdd_p = SAFEDIV_POS(nm_px, -r->ctx_tau);
  }
  if (r->bty_tau < -tol) {
    nm_aty = NORM(r->aty, n);
    r->res_infeas = SAFEDIV_POS(nm_aty, -r->bty_tau);
  }
}

static void unnormalize_residuals(ScsWork *w) {
  ScsResiduals *r_n = w->r_normalized; /* normalized residuals */
  ScsResiduals *r = w->r_orig;         /* original problem residuals */
  scs_float pd = w->scal->primal_scale * w->scal->dual_scale;

  /* copy vars */
  r->last_iter = r_n->last_iter;
  r->tau = r_n->tau;

  /* unnormalize */
  r->kap = r_n->kap / pd;
  r->bty_tau = r_n->bty_tau / pd;
  r->ctx_tau = r_n->ctx_tau / pd;
  r->xt_p_x_tau = r_n->xt_p_x_tau / pd;
  r->xt_p_x = r_n->xt_p_x / pd;
  r->ctx = r_n->ctx / pd;
  r->bty = r_n->bty / pd;
  r->pobj = r_n->pobj / pd;
  r->dobj = r_n->dobj / pd;
  r->gap = r_n->gap / pd;

  /* Fuse the six memcpy+un_normalize calls into two loops.
   * Primal: divide by D[i]*dual_scale. Dual: divide by E[i]*primal_scale.
   * This reduces 6 memcpy + 6 scale passes to 2 passes. */
  {
    scs_int i;
    const scs_float *D = w->scal->D, *E = w->scal->E;
    scs_float inv_ds = 1.0 / (w->scal->dual_scale);
    scs_float inv_ps = 1.0 / (w->scal->primal_scale);
    for (i = 0; i < w->d->m; ++i) {
      scs_float f = inv_ds / D[i];
      r->ax[i]        = r_n->ax[i]        * f;
      r->ax_s[i]      = r_n->ax_s[i]      * f;
      r->ax_s_btau[i] = r_n->ax_s_btau[i] * f;
    }
    for (i = 0; i < w->d->n; ++i) {
      scs_float f = inv_ps / E[i];
      r->aty[i]         = r_n->aty[i]         * f;
      r->px[i]          = r_n->px[i]           * f;
      r->px_aty_ctau[i] = r_n->px_aty_ctau[i]  * f;
    }
  }

  compute_residuals(r, w->d->m, w->d->n, pd);
}

/* calculates un-normalized residual quantities */
/* this is somewhat slow but not a bottleneck */
static void populate_residual_struct(ScsWork *w, scs_int iter) {
  scs_int i, n = w->d->n, m = w->d->m;
  /* normalized x,y,s terms */
  scs_float *x = w->xys_normalized->x;
  scs_float *y = w->xys_normalized->y;
  scs_float *s = w->xys_normalized->s;
  ScsResiduals *r = w->r_normalized; /* normalized residuals */

  /* checks if the residuals are unchanged by checking iteration */
  if (r->last_iter == iter) {
    return;
  }
  r->last_iter = iter;

  memcpy(x, w->u, n * sizeof(scs_float));
  memcpy(y, &(w->u[n]), m * sizeof(scs_float));
  memcpy(s, &(w->rsk[n]), m * sizeof(scs_float));

  r->tau = ABS(w->u[n + m]);
  r->kap = ABS(w->rsk[n + m]);

  /**************** PRIMAL *********************/
  memset(r->ax, 0, m * sizeof(scs_float));
  /* ax = Ax */
  SCS(accum_by_a)(w->d->A, x, r->ax);

  /* Build ax_s and ax_s_btau in one fused pass (saves 2 memcpy + 2 axpy). */
  for (i = 0; i < m; ++i) {
    r->ax_s[i] = r->ax[i] + s[i];
    r->ax_s_btau[i] = r->ax_s[i] - r->tau * w->d->b[i];
  }

  /**************** DUAL *********************/
  memset(r->px, 0, n * sizeof(scs_float));
  if (w->d->P) {
    /* px = Px */
    SCS(accum_by_p)(w->d->P, x, r->px);
    r->xt_p_x_tau = SCS(dot)(r->px, x, n);
  } else {
    r->xt_p_x_tau = 0.;
  }

  memset(r->aty, 0, n * sizeof(scs_float));
  /* aty = A'y */
  SCS(accum_by_atrans)(w->d->A, y, r->aty);

  /* Build px_aty_ctau in one fused pass (saves 1 memcpy + 2 axpy). */
  for (i = 0; i < n; ++i) {
    r->px_aty_ctau[i] = r->px[i] + r->aty[i] + r->tau * w->d->c[i];
  }

  /**************** OTHERS *****************/
  r->bty_tau = SCS(dot)(y, w->d->b, m);
  r->ctx_tau = SCS(dot)(x, w->d->c, n);

  r->bty = SAFEDIV_POS(r->bty_tau, r->tau);
  r->ctx = SAFEDIV_POS(r->ctx_tau, r->tau);
  r->xt_p_x = SAFEDIV_POS(r->xt_p_x_tau, r->tau * r->tau);

  r->gap = ABS(r->xt_p_x + r->ctx + r->bty);
  r->pobj = r->xt_p_x / 2. + r->ctx;
  r->dobj = -r->xt_p_x / 2. - r->bty;

  compute_residuals(r, m, n, 1.0);

  if (w->stgs->normalize) {
    memcpy(w->xys_orig->x, w->xys_normalized->x, n * sizeof(scs_float));
    memcpy(w->xys_orig->y, w->xys_normalized->y, m * sizeof(scs_float));
    memcpy(w->xys_orig->s, w->xys_normalized->s, m * sizeof(scs_float));
    SCS(un_normalize_sol)(w->scal, w->xys_orig);
    unnormalize_residuals(w);
  }
}

/* ==================== Convergence Checking ========================= */

static scs_int has_converged(ScsWork *w, scs_int iter) {
  scs_float abs_xt_p_x, abs_ctx, abs_bty;
  scs_float nm_s, nm_px, nm_aty, nm_ax;
  scs_float grl, prl, drl;
  scs_float eps_abs = w->stgs->eps_abs;
  scs_float eps_rel = w->stgs->eps_rel;
  scs_float eps_infeas = w->stgs->eps_infeas;

  ScsResiduals *r = w->r_orig;

  if (r->tau > 0.) {
    abs_xt_p_x = ABS(r->xt_p_x);
    abs_ctx = ABS(r->ctx);
    abs_bty = ABS(r->bty);

    nm_s = NORM(w->xys_orig->s, w->d->m);
    nm_px = NORM(r->px, w->d->n);
    nm_aty = NORM(r->aty, w->d->n);
    nm_ax = NORM(r->ax, w->d->m);
    /* xt_p_x, ctx, bty already have tau divided out */
    grl = MAX(MAX(abs_xt_p_x, abs_ctx), abs_bty);
    /* s, ax, px, aty do *not* have tau divided out, so need to divide */
    prl = MAX(MAX(w->nm_b_orig * r->tau, nm_s), nm_ax) / r->tau;
    drl = MAX(MAX(w->nm_c_orig * r->tau, nm_px), nm_aty) / r->tau;
    if (isless(r->res_pri, eps_abs + eps_rel * prl) &&
        isless(r->res_dual, eps_abs + eps_rel * drl) &&
        isless(r->gap, eps_abs + eps_rel * grl)) {
      return SCS_SOLVED;
    }
  }
  if (isless(r->res_unbdd_a, eps_infeas) &&
      isless(r->res_unbdd_p, eps_infeas)) {
    return SCS_UNBOUNDED;
  }
  if (isless(r->res_infeas, eps_infeas)) {
    return SCS_INFEASIBLE;
  }
  return 0;
}

/* =================== Warm / Cold Start Helpers ==================== */

static inline scs_int _is_nan(scs_float x) {
  return x != x;
}

/* given x,y,s warm start, set v = [x; s / R + y; 1]
 * check for nans and set to zero if present
 */
static void warm_start_vars(ScsWork *w, ScsSolution *sol) {
  scs_int n = w->d->n, m = w->d->m, i;
  scs_float *v = w->v;
  /* normalize the warm-start */
  if (w->stgs->normalize) {
    SCS(normalize_sol)(w->scal, sol);
  }
  for (i = 0; i < n; ++i) {
    v[i] = _is_nan(sol->x[i]) ? 0. : sol->x[i];
  }
  for (i = 0; i < m; ++i) {
    v[i + n] = sol->y[i] + sol->s[i] / w->diag_r[i + n];
    v[i + n] = _is_nan(v[i + n]) ? 0. : v[i + n];
  }
  v[n + m] = 1.0; /* tau = 1 */
  /* un-normalize so sol unchanged */
  if (w->stgs->normalize) {
    SCS(un_normalize_sol)(w->scal, sol);
  }
}

static void cold_start_vars(ScsWork *w) {
  scs_int l = w->d->n + w->d->m + 1;
  memset(w->v, 0, l * sizeof(scs_float));
  w->v[l - 1] = 1.;
}

/* ====================== ADMM Iteration ============================= */

static scs_float root_plus(ScsWork *w, scs_float *p, scs_float *mu,
                           scs_float eta) {
  /* Compute all five weighted dot products (g'Rg, mu'Rg, p'Rg, p'Rp, p'Rmu)
   * in a single pass over diag_r to minimise memory traffic. */
  scs_int i, nm = w->d->n + w->d->m;
  scs_float gg = 0., mug = 0., pg = 0., pp = 0., pmu = 0.;
  scs_float a, b, c, rad, tau_scale = w->diag_r[nm];
  const scs_float *g = w->g, *r = w->diag_r;
  for (i = 0; i < nm; ++i) {
    scs_float ri = r[i], gi = g[i], pi = p[i], mui = mu[i];
    gg  += gi  * gi  * ri;
    mug += mui * gi  * ri;
    pg  += pi  * gi  * ri;
    pp  += pi  * pi  * ri;
    pmu += pi  * mui * ri;
  }
  a = tau_scale + gg;
  b = mug - 2 * pg - eta * tau_scale;
  c = pp - pmu;
  rad = b * b - 4 * a * c;
  return (-b + SQRTF(MAX(rad, 0.))) / (2 * a);
}

/* status != 0 indicates failure */
static scs_int project_lin_sys(ScsWork *w, scs_int iter) {
  scs_int n = w->d->n, m = w->d->m, l = n + m + 1, status, i;
  scs_float *warm_start = SCS_NULL;
  scs_float tol = -1.0; /* only used for indirect methods, overridden later */
  /* Copy and scale in one pass, eliminating the intermediate memcpy. */
  for (i = 0; i < n; ++i) {
    w->u_t[i] = w->v[i] * w->diag_r[i];
  }
  for (i = n; i < l - 1; ++i) {
    w->u_t[i] = -w->v[i] * w->diag_r[i];
  }
  w->u_t[l - 1] = w->v[l - 1];
#if INDIRECT > 0
  scs_float nm_ax_s_btau, nm_px_aty_ctau, nm_ws;
  /* compute warm start using the cone projection output */
  nm_ax_s_btau = CG_NORM(w->r_normalized->ax_s_btau, m);
  nm_px_aty_ctau = CG_NORM(w->r_normalized->px_aty_ctau, n);
  warm_start = w->lin_sys_warm_start;
  /* warm_start = u[:n] + tau * g[:n] */
  memcpy(warm_start, w->u, n * sizeof(scs_float));
  SCS(add_scaled_array)(warm_start, w->g, n, w->u[l - 1]);
  /* use normalized residuals to compute tolerance */
  tol = MIN(nm_ax_s_btau, nm_px_aty_ctau);
  /* tol ~ O(1/k^(1+eps)) guarantees convergence */
  /* use warm-start to calculate tolerance rather than w->u_t, since warm_start
   * should be approximately equal to the true solution */
  nm_ws = CG_NORM(warm_start, n) / POWF((scs_float)iter + 1, CG_RATE);
  tol = CG_TOL_FACTOR * MIN(tol, nm_ws);
  tol = MAX(CG_BEST_TOL, tol);
#endif
  status = scs_solve_lin_sys(w->p, w->u_t, warm_start, tol);
  if (iter < FEASIBLE_ITERS) {
    w->u_t[l - 1] = 1.;
  } else {
    w->u_t[l - 1] = root_plus(w, w->u_t, w->v, w->v[l - 1]);
  }
  SCS(add_scaled_array)(w->u_t, w->g, l - 1, -w->u_t[l - 1]);
  return status;
}

/* Compute the [r;s;kappa] iterate
 *
 *  rsk^{k+1} = R ( u^{k+1} + v^k - 2 * u_t^{k+1} )
 *
 *  uses Moreau decomposition to get projection onto dual cone
 *  since it depends on v^k MUST be called before update_dual_vars is done
 *  (no effect of w->stgs->alpha here).
 */
static void compute_rsk(ScsWork *w) {
  scs_int i, l = w->d->m + w->d->n + 1;
  for (i = 0; i < l; ++i) {
    w->rsk[i] = (w->v[i] + w->u[i] - 2 * w->u_t[i]) * w->diag_r[i];
  }
}

static void update_dual_vars(ScsWork *w) {
  scs_int i, l = w->d->n + w->d->m + 1;
  for (i = 0; i < l; ++i) {
    w->v[i] += w->stgs->alpha * (w->u[i] - w->u_t[i]);
  }
}

/* status < 0 indicates failure */
static scs_int project_cones(ScsWork *w, const ScsCone *k, scs_int iter) {
  scs_int i, n = w->d->n, l = w->d->n + w->d->m + 1, status;
  for (i = 0; i < l; ++i) {
    w->u[i] = 2 * w->u_t[i] - w->v[i];
  }
  /* u = [x;y;tau] */
  status =
      SCS(proj_dual_cone)(&(w->u[n]), w->cone_work, w->scal, &(w->diag_r[n]));
  if (iter < FEASIBLE_ITERS) {
    w->u[l - 1] = 1.0;
  } else {
    w->u[l - 1] = MAX(w->u[l - 1], 0.);
  }
  return status;
}

/* scs is homogeneous so scale the iterate to keep norm reasonable */
static inline void normalize_v(scs_float *v, scs_int len) {
  scs_float v_norm = SCS(norm_2)(v, len); /* always l2 norm */
  if (v_norm == 0.) {
    scs_printf("WARNING: normalize_v called with zero-norm iterate; this is "
               "highly pathological (e.g., strong duality may not hold).\n");
    return;
  }
  SCS(scale_array)(v, SQRTF((scs_float)len) * ITERATE_NORM / v_norm, len);
}

/* ================== Solution Extraction / Finalization ============== */

static void sety(const ScsWork *w, ScsSolution *sol) {
  if (!sol->y) {
    sol->y = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  }
  memcpy(sol->y, &(w->u[w->d->n]), w->d->m * sizeof(scs_float));
}

/* s is contained in rsk */
static void sets(const ScsWork *w, ScsSolution *sol) {
  if (!sol->s) {
    sol->s = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  }
  memcpy(sol->s, &(w->rsk[w->d->n]), w->d->m * sizeof(scs_float));
}

static void setx(const ScsWork *w, ScsSolution *sol) {
  if (!sol->x) {
    sol->x = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
  }
  memcpy(sol->x, w->u, w->d->n * sizeof(scs_float));
}

static void set_solved(const ScsWork *w, ScsSolution *sol, ScsInfo *info) {
  SCS(scale_array)(sol->x, SAFEDIV_POS(1.0, w->r_orig->tau), w->d->n);
  SCS(scale_array)(sol->y, SAFEDIV_POS(1.0, w->r_orig->tau), w->d->m);
  SCS(scale_array)(sol->s, SAFEDIV_POS(1.0, w->r_orig->tau), w->d->m);
  info->gap = w->r_orig->gap;
  info->res_pri = w->r_orig->res_pri;
  info->res_dual = w->r_orig->res_dual;
  info->pobj = w->r_orig->xt_p_x / 2. + w->r_orig->ctx;
  info->dobj = -w->r_orig->xt_p_x / 2. - w->r_orig->bty;
  strcpy(info->status, "solved");
  info->status_val = SCS_SOLVED;
}

static void set_infeasible(const ScsWork *w, ScsSolution *sol, ScsInfo *info) {
  SCS(scale_array)(sol->y, -1 / w->r_orig->bty_tau, w->d->m);
  SCS(scale_array)(sol->x, NAN, w->d->n);
  SCS(scale_array)(sol->s, NAN, w->d->m);
  info->gap = NAN;
  info->res_pri = NAN;
  info->res_dual = NAN;
  info->pobj = INFINITY;
  info->dobj = INFINITY;
  strcpy(info->status, "infeasible");
  info->status_val = SCS_INFEASIBLE;
}

static void set_unbounded(const ScsWork *w, ScsSolution *sol, ScsInfo *info) {
  SCS(scale_array)(sol->x, -1 / w->r_orig->ctx_tau, w->d->n);
  SCS(scale_array)(sol->s, -1 / w->r_orig->ctx_tau, w->d->m);
  SCS(scale_array)(sol->y, NAN, w->d->m);
  info->gap = NAN;
  info->res_pri = NAN;
  info->res_dual = NAN;
  info->pobj = -INFINITY;
  info->dobj = -INFINITY;
  strcpy(info->status, "unbounded");
  info->status_val = SCS_UNBOUNDED;
}

/* not yet converged, take best guess */
static void set_unfinished(const ScsWork *w, ScsSolution *sol, ScsInfo *info) {
  if (w->r_orig->kap > w->r_orig->tau &&
      (w->r_orig->bty_tau < 0 || w->r_orig->ctx_tau < 0)) {
    if (w->r_orig->bty_tau < 0 &&
        w->r_orig->bty_tau < w->r_orig->ctx_tau) {
      set_infeasible(w, sol, info);
      info->status_val = SCS_INFEASIBLE_INACCURATE;
    } else {
      set_unbounded(w, sol, info);
      info->status_val = SCS_UNBOUNDED_INACCURATE;
    }
  } else if (w->r_orig->tau > 0) {
    set_solved(w, sol, info);
    info->status_val = SCS_SOLVED_INACCURATE;
  } else {
    scs_printf("ERROR: could not determine problem status.\n");
    info->status_val = SCS_FAILED;
  }
  /* Append inaccurate to the status string */
  if (w->time_limit_reached) {
    strcat(info->status, " (inaccurate - reached time_limit_secs)");
  } else if (info->iter >= w->stgs->max_iters) {
    strcat(info->status, " (inaccurate - reached max_iters)");
  } else {
    scs_printf("ERROR: should not be in this state (1).\n");
  }
}

/* sets solutions, re-scales by inner prods if infeasible or unbounded */
static void finalize(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                     scs_int iter) {
  scs_float nm_s, nm_y, sty;
  setx(w, sol);
  sety(w, sol);
  sets(w, sol);
  if (w->stgs->normalize) {
    SCS(un_normalize_sol)(w->scal, sol);
  }
  populate_residual_struct(w, iter);

  nm_s = SCS(norm_inf)(sol->s, w->d->m);
  nm_y = SCS(norm_inf)(sol->y, w->d->m);
  sty = SCS(dot)(sol->s, sol->y, w->d->m);

  info->setup_time = w->setup_time;
  info->iter = iter;
  info->res_infeas = w->r_orig->res_infeas;
  info->res_unbdd_a = w->r_orig->res_unbdd_a;
  info->res_unbdd_p = w->r_orig->res_unbdd_p;
  info->scale = w->stgs->scale;
  info->scale_updates = w->scale_updates;
  info->rejected_accel_steps = w->rejected_accel_steps;
  info->accepted_accel_steps = w->accepted_accel_steps;
  info->comp_slack = ABS(sty);
#ifdef SPECTRAL_TIMING_FLAG
  info->ave_time_matrix_cone_proj = w->cone_work->tot_time_mat_cone_proj / iter;
  info->ave_time_vector_cone_proj = w->cone_work->tot_time_vec_cone_proj / iter;
#endif
  if (info->comp_slack > 1e-5 * MAX(nm_s, nm_y)) {
    scs_printf("WARNING - large complementary slackness residual: %f\n",
               info->comp_slack);
  }
  switch (info->status_val) {
  case SCS_SOLVED:
    set_solved(w, sol, info);
    break;
  case SCS_INFEASIBLE:
    set_infeasible(w, sol, info);
    break;
  case SCS_UNBOUNDED:
    set_unbounded(w, sol, info);
    break;
  case SCS_UNFINISHED: /* When SCS reaches max_iters or time_limit_secs */
    set_unfinished(w, sol, info);
    break;
  default:
    scs_printf("ERROR: should not be in this state (2).\n");
  }
}

/* ================ Workspace Init / Scale Updating ================== */

/* Sets the diag_r vector, given the scale parameters in work */
static void set_diag_r(ScsWork *w) {
  scs_int i;
  for (i = 0; i < w->d->n; ++i) {
    w->diag_r[i] = w->stgs->rho_x;
  }
  /* use cone information to set R_y */
  SCS(set_r_y)(w->cone_work, w->stgs->scale, &(w->diag_r[w->d->n]));
  /* if modified need to SCS(enforce_cone_boundaries)(...) */
  w->diag_r[w->d->n + w->d->m] = TAU_FACTOR;
}

static ScsWork *init_work(const ScsData *d, const ScsCone *k,
                          const ScsSettings *stgs) {
  ScsWork *w = (ScsWork *)scs_calloc(1, sizeof(ScsWork));
  scs_int l = d->n + d->m + 1;
  if (stgs->verbose) {
    print_init_header(d, k, stgs);
  }
  if (!w) {
    scs_printf("ERROR: allocating work failure\n");
    return SCS_NULL;
  }
  /* deep copy data */
  w->d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  if (!w->d || !SCS(deep_copy_data)(w->d, d)) {
    scs_printf("ERROR: data copy failure\n");
    scs_finish(w);
    return SCS_NULL;
  }
  d = SCS_NULL; /* for safety */

  /* deep copy cone */
  w->k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  if (!w->k || !SCS(deep_copy_cone)(w->k, k)) {
    scs_printf("ERROR: cone copy failure\n");
    scs_finish(w);
    return SCS_NULL;
  }
  k = SCS_NULL; /* for safety */

  /* deep copy settings */
  w->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  if (!w->stgs || !SCS(deep_copy_stgs)(w->stgs, stgs)) {
    scs_printf("ERROR: settings copy failure\n");
    scs_finish(w);
    return SCS_NULL;
  }
  stgs = SCS_NULL; /* for safety */

  /* allocate workspace: */
  w->u = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->u_t = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v_prev = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->rsk = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->h = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->g = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->lin_sys_warm_start = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
  w->diag_r = (scs_float *)scs_calloc(l, sizeof(scs_float));
  /* x,y,s struct */
  w->xys_orig = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  w->xys_orig->x = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
  w->xys_orig->s = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  w->xys_orig->y = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  w->r_orig = init_residuals(w->d);
  w->b_orig = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  w->c_orig = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));

  if (!w->u || !w->u_t || !w->v || !w->v_prev || !w->rsk || !w->h || !w->g ||
      !w->lin_sys_warm_start || !w->diag_r || !w->xys_orig ||
      !w->xys_orig->x || !w->xys_orig->s || !w->xys_orig->y || !w->r_orig ||
      !w->b_orig || !w->c_orig) {
    scs_printf("ERROR: work memory allocation failure\n");
    scs_finish(w);
    return SCS_NULL;
  }

  if (!(w->cone_work = SCS(init_cone)(w->k, w->d->m))) {
    scs_printf("ERROR: init_cone failure\n");
    scs_finish(w);
    return SCS_NULL;
  }
  set_diag_r(w);

  if (w->stgs->normalize) {
    w->xys_normalized = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
    w->xys_normalized->x = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
    w->xys_normalized->s = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
    w->xys_normalized->y = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
    w->r_normalized = init_residuals(w->d);
    if (!w->xys_normalized || !w->xys_normalized->x || !w->xys_normalized->s ||
        !w->xys_normalized->y || !w->r_normalized) {
      scs_printf("ERROR: normalized work memory allocation failure\n");
      scs_finish(w);
      return SCS_NULL;
    }
    /* this allocates memory that must be freed */
    w->scal = SCS(normalize_a_p)(w->d->P, w->d->A, w->cone_work);
    if (!w->scal) {
      scs_printf("ERROR: normalize_a_p failure\n");
      scs_finish(w);
      return SCS_NULL;
    }
  } else {
    w->xys_normalized = w->xys_orig;
    w->r_normalized = w->r_orig;
    w->scal = SCS_NULL;
  }
  /* set w->*_orig and performs normalization if appropriate */
  scs_update(w, w->d->b, w->d->c);

  if (!(w->p = scs_init_lin_sys_work(w->d->A, w->d->P, w->diag_r))) {
    scs_printf("ERROR: init_lin_sys_work failure\n");
    scs_finish(w);
    return SCS_NULL;
  }
  if (w->stgs->acceleration_lookback) {
    /* TODO(HACK!) negative acceleration_lookback interpreted as type-II */
    if (!(w->accel = aa_init(l, IABS(w->stgs->acceleration_lookback),
                             w->stgs->acceleration_lookback > 0,
                             w->stgs->acceleration_lookback > 0
                                 ? AA_REGULARIZATION_TYPE_1
                                 : AA_REGULARIZATION_TYPE_2,
                             AA_RELAXATION, AA_SAFEGUARD_FACTOR,
                             AA_MAX_WEIGHT_NORM, VERBOSITY))) {
      if (w->stgs->verbose) {
        scs_printf("WARN: aa_init returned NULL, no acceleration applied.\n");
      }
    }
  } else {
    w->accel = SCS_NULL;
  }
  return w;
}

static void update_work_cache(ScsWork *w) {
  /* g = (I + M)^{-1} [c; -b]
   * Build g = [c; -b] directly in one pass, avoiding a separate memcpy of h
   * followed by a negate pass over g[n:]. */
  scs_int i, n = w->d->n, m = w->d->m;
  memcpy(w->g, w->d->c, n * sizeof(scs_float));
  for (i = 0; i < m; ++i) {
    w->g[n + i] = -w->d->b[i];
  }
  scs_solve_lin_sys(w->p, w->g, SCS_NULL, CG_BEST_TOL);
}

/* Reset quantities specific to current solve */
static void reset_tracking(ScsWork *w) {
  w->last_scale_update_iter = 0;
  w->sum_log_scale_factor = 0.;
  w->n_log_scale_factor = 0;
  w->scale_updates = 0;
  w->time_limit_reached = 0;
  /* Acceleration */
  w->rejected_accel_steps = 0;
  w->accepted_accel_steps = 0;
  w->aa_norm = 0.;
  /* Need this to force residual calc if previous solve solved at iter 0 */
  w->r_normalized->last_iter = -1;
  w->r_orig->last_iter = -1;
}

static scs_int update_work(ScsWork *w, ScsSolution *sol) {
  reset_tracking(w);

  if (w->stgs->warm_start) {
    warm_start_vars(w, sol);
  } else {
    cold_start_vars(w);
  }

  update_work_cache(w);
  return 0;
}

/* will update if the factor is outside of range */
static scs_int should_update_r(scs_float factor) {
  return (factor > SQRTF(10.) || factor < 1. / SQRTF(10.));
}

static scs_int update_scale(ScsWork *w, const ScsCone *k, scs_int iter) {
  scs_int i;
  scs_float factor, new_scale, relative_res_pri, relative_res_dual;
  scs_float denom_pri, denom_dual;

  ScsResiduals *r = w->r_orig;

  scs_float nm_ax = SCALE_NORM(r->ax, w->d->m);
  scs_float nm_s = SCALE_NORM(w->xys_orig->s, w->d->m);
  scs_float nm_px_aty_ctau = SCALE_NORM(r->px_aty_ctau, w->d->n);
  scs_float nm_px = SCALE_NORM(r->px, w->d->n);
  scs_float nm_aty = SCALE_NORM(r->aty, w->d->n);
  scs_float nm_ax_s_btau = SCALE_NORM(r->ax_s_btau, w->d->m);

  scs_int iters_since_last_update = iter - w->last_scale_update_iter;
  /* ||Ax + s - b * tau|| */
  denom_pri = MAX(nm_ax, nm_s);
  denom_pri = MAX(denom_pri, w->nm_b_orig * r->tau);
  relative_res_pri = SAFEDIV_POS(nm_ax_s_btau, denom_pri);
  /* ||Px + A'y + c * tau|| */
  denom_dual = MAX(nm_px, nm_aty);
  denom_dual = MAX(denom_dual, w->nm_c_orig * r->tau);
  relative_res_dual = SAFEDIV_POS(nm_px_aty_ctau, denom_dual);

  /* higher scale makes res_pri go down faster, so increase if res_pri larger */
  /* clamp to avoid log(0) which would NaN-poison sum_log_scale_factor */
  relative_res_pri = MAX(relative_res_pri, _DIV_EPS_TOL);
  relative_res_dual = MAX(relative_res_dual, _DIV_EPS_TOL);
  w->sum_log_scale_factor += log(relative_res_pri) - log(relative_res_dual);
  w->n_log_scale_factor++;

  /* geometric mean */
  factor =
      SQRTF(exp(w->sum_log_scale_factor / (scs_float)(w->n_log_scale_factor)));

  /* need at least RESCALING_MIN_ITERS since last update */
  if (iters_since_last_update < RESCALING_MIN_ITERS) {
    return 0;
  }
  new_scale =
      MIN(MAX(w->stgs->scale * factor, MIN_SCALE_VALUE), MAX_SCALE_VALUE);
  if (new_scale == w->stgs->scale) {
    return 0;
  }
  if (should_update_r(factor)) {
    scs_int linsys_status;
    w->scale_updates++;
    w->sum_log_scale_factor = 0;
    w->n_log_scale_factor = 0;
    w->last_scale_update_iter = iter;
    w->stgs->scale = new_scale;

    /* update diag r vector */
    set_diag_r(w);

    /* update linear systems */
    linsys_status = scs_update_lin_sys_diag_r(w->p, w->diag_r);
    if (linsys_status < 0) {
      return linsys_status;
    }

    /* update pre-solved quantities */
    update_work_cache(w);

    /* reset acceleration so that old iterates aren't affecting new values */
    if (w->accel) {
      aa_reset(w->accel);
    }
    /* update v, using fact that rsk, u, u_t vectors should be the same */
    /* solve: R^+ (v^+ + u - 2u_t) = rsk = R(v + u - 2u_t)
     *  => v^+ = R+^-1 rsk + 2u_t - u
     */
    for (i = 0; i < w->d->n + w->d->m + 1; i++) {
      w->v[i] = w->rsk[i] / w->diag_r[i] + 2 * w->u_t[i] - w->u[i];
    }
  }
  return 0;
}

/* ========================== Public API ============================= */

ScsWork *scs_init(const ScsData *d, const ScsCone *k, const ScsSettings *stgs) {
  ScsWork *w;
  SCS(timer) init_timer;
  if (!d || !k || !stgs) {
    scs_printf("ERROR: Missing ScsData, ScsCone, or ScsSettings input\n");
    return SCS_NULL;
  }
#ifdef SCS_MKL
  if (scs_init_mkl_runtime() != 0) {
    scs_printf("ERROR: MKL runtime initialization failed before solver setup. "
               "See the MKL message above for the specific cause.\n");
    return SCS_NULL;
  }
#endif
#if NO_VALIDATE == 0
  if (validate(d, k, stgs) < 0) {
    scs_printf("ERROR: Validation returned failure\n");
    return SCS_NULL;
  }
#endif
  scs_start_interrupt_listener();
#if VERBOSITY > 0
  scs_printf("size of scs_int = %lu, size of scs_float = %lu\n",
             (unsigned long)sizeof(scs_int), (unsigned long)sizeof(scs_float));
#endif
  SCS(tic)(&init_timer);
  if (stgs->write_data_filename) {
    scs_printf("Writing raw problem data to %s\n", stgs->write_data_filename);
    SCS(write_data)(d, k, stgs);
  }
  if (stgs->log_csv_filename) {
    scs_printf("Logging run data to %s\n", stgs->log_csv_filename);
    /* logging done every iteration */
  }
  w = init_work(d, k, stgs);
  if (w) {
    w->setup_time = SCS(tocq)(&init_timer);
  }
  scs_end_interrupt_listener();
  return w;
}

scs_int scs_update(ScsWork *w, scs_float *b, scs_float *c) {
  SCS(timer) update_timer;
  SCS(tic)(&update_timer);

  if (b) {
    if (w->b_orig != b) {
      memcpy(w->b_orig, b, w->d->m * sizeof(scs_float));
    }
    if (w->d->b != b) {
      memcpy(w->d->b, b, w->d->m * sizeof(scs_float));
    }
    w->nm_b_orig = NORM(w->b_orig, w->d->m);
  } else {
    /* b_orig unchanged so no need to recompute the norm */
    memcpy(w->d->b, w->b_orig, w->d->m * sizeof(scs_float));
  }

  if (c) {
    if (w->c_orig != c) {
      memcpy(w->c_orig, c, w->d->n * sizeof(scs_float));
    }
    if (w->d->c != c) {
      memcpy(w->d->c, c, w->d->n * sizeof(scs_float));
    }
    w->nm_c_orig = NORM(w->c_orig, w->d->n);
  } else {
    /* c_orig unchanged so no need to recompute the norm */
    memcpy(w->d->c, w->c_orig, w->d->n * sizeof(scs_float));
  }

  /* normalize */
  if (w->scal) {
    SCS(normalize_b_c)(w->scal, w->d->b, w->d->c);
  }

  /* override setup time with update time, since the update is the 'setup' */
  w->setup_time = SCS(tocq)(&update_timer);
  return 0;
}

scs_int scs_solve(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                  scs_int warm_start) {
  scs_int i;
  SCS(timer) solve_timer, lin_sys_timer, cone_timer, accel_timer;
  scs_float total_accel_time = 0.0, total_cone_time = 0.0,
            total_lin_sys_time = 0.0;
  if (!sol || !w || !info) {
    scs_printf("ERROR: missing ScsWork, ScsSolution or ScsInfo input\n");
    return SCS_FAILED;
  }
  scs_int l = w->d->m + w->d->n + 1;
  const ScsCone *k = w->k;
  ScsSettings *stgs = w->stgs;
  /* set warm start */
  stgs->warm_start = warm_start;

  /* initialize ctrl-c support */
  scs_start_interrupt_listener();
  SCS(tic)(&solve_timer);
  strcpy(info->lin_sys_solver, scs_get_lin_sys_method());
  info->status_val = SCS_UNFINISHED; /* not yet converged */
  update_work(w, sol);

  if (w->stgs->verbose) {
    print_header(w, k);
  }

  /* SCS */
  for (i = 0; i < w->stgs->max_iters; ++i) {
    /* Accelerate here so that last step always projection onto cone */
    /* this ensures the returned iterates always satisfy conic constraints */
    if (w->accel) {
      SCS(tic)(&accel_timer);
      if (i > 0 && i % w->stgs->acceleration_interval == 0) {
        /* v overwritten with AA output here */
        w->aa_norm = aa_apply(w->v, w->v_prev, w->accel);
      }
      total_accel_time += SCS(tocq)(&accel_timer);
    }

    if (i >= FEASIBLE_ITERS) {
      /* normalize v *after* applying any acceleration */
      /* the input to the DR step should be normalized */
      normalize_v(w->v, l);
    }

    /* store v_prev = v for AA safeguard; skip when acceleration is off */
    if (w->accel) {
      memcpy(w->v_prev, w->v, l * sizeof(scs_float));
    }

    /******************* linear system solve ********************/
    SCS(tic)(&lin_sys_timer);
    if (project_lin_sys(w, i) != 0) {
      return failure(w, w->d->m, w->d->n, sol, info, SCS_FAILED,
                     "error in project_lin_sys", "failure");
    }
    total_lin_sys_time += SCS(tocq)(&lin_sys_timer);

    /****************** project onto the cones ******************/
    SCS(tic)(&cone_timer);
    if (project_cones(w, k, i) < 0) {
      return failure(w, w->d->m, w->d->n, sol, info, SCS_FAILED,
                     "error in project_cones", "failure");
    }
    total_cone_time += SCS(tocq)(&cone_timer);

    /* compute [r;s;kappa], must be before dual var update */
    /* since Moreau decomp logic relies on v at start */
    compute_rsk(w);

    if (i % CONVERGED_INTERVAL == 0) {
      if (scs_is_interrupted()) {
        return failure(w, w->d->m, w->d->n, sol, info, SCS_SIGINT,
                       "interrupted", "interrupted");
      }
      populate_residual_struct(w, i);
      if ((info->status_val = has_converged(w, i)) != 0) {
        break;
      }
      if (w->stgs->time_limit_secs) {
        if (SCS(tocq)(&solve_timer) > 1000. * w->stgs->time_limit_secs) {
          w->time_limit_reached = 1;
          break;
        }
      }
    }

    /* Compute residuals. */
    if (w->stgs->verbose && i % PRINT_INTERVAL == 0) {
      populate_residual_struct(w, i);
      print_summary(w, i, &solve_timer);
    }

    /* If residuals are fresh then maybe compute new scale. */
    if (w->stgs->adaptive_scale && i == w->r_orig->last_iter) {
      if (update_scale(w, k, i) < 0) {
        return failure(w, w->d->m, w->d->n, sol, info, SCS_FAILED,
                       "error in update_scale", "failure");
      }
    }

    /****************** dual variable step **********************/
    /* do this after update_scale due to remapping that happens there */
    update_dual_vars(w);

    /* AA safeguard check.
     * Perform safeguarding *after* convergence check to prevent safeguard
     * overwriting converged iterate, since safeguard is on `v` and convergence
     * is on `u`.
     */
    if (w->accel && i % w->stgs->acceleration_interval == 0 && w->aa_norm > 0) {
      if (aa_safeguard(w->v, w->v_prev, w->accel) < 0) {
        /* TODO should we copy u from u_prev here too? Then move above, possibly
         * better residual calculation and scale updating. */
        w->rejected_accel_steps++;
      } else {
        w->accepted_accel_steps++;
      }
    }

    /* Log *after* updating scale so residual recalc does not affect alg */
    if (w->stgs->log_csv_filename) {
      /* calc residuals every iter if logging to csv */
      populate_residual_struct(w, i);
      SCS(log_data_to_csv)(k, stgs, w, i, &solve_timer);
    }
  }

  /* Final logging after full run */
  if (w->stgs->log_csv_filename) {
    populate_residual_struct(w, i);
    SCS(log_data_to_csv)(k, stgs, w, i, &solve_timer);
  }

  if (w->stgs->verbose) {
    populate_residual_struct(w, i);
    print_summary(w, i, &solve_timer);
  }

  /* populate solution vectors (unnormalized) and info */
  finalize(w, sol, info, i);

  /* populate timings */
  info->solve_time = SCS(tocq)(&solve_timer);
  info->lin_sys_time = total_lin_sys_time;
  info->cone_time = total_cone_time;
  info->accel_time = total_accel_time;

  if (w->stgs->verbose) {
    print_footer(info);
  }

  scs_end_interrupt_listener();
  return info->status_val;
}

void scs_finish(ScsWork *w) {
  if (w) {
    if (w->cone_work) {
      SCS(finish_cone)(w->cone_work);
    }
    if (w->p) {
      scs_free_lin_sys_work(w->p);
    }
    if (w->accel) {
      aa_finish(w->accel);
    }
    scs_free(w->u);
    scs_free(w->u_t);
    scs_free(w->v);
    scs_free(w->v_prev);
    scs_free(w->rsk);
    scs_free(w->h);
    scs_free(w->g);
    scs_free(w->b_orig);
    scs_free(w->c_orig);
    scs_free(w->lin_sys_warm_start);
    scs_free(w->diag_r);
    SCS(free_sol)(w->xys_orig);
    if (w->scal) {
      scs_free(w->scal->D);
      scs_free(w->scal->E);
      scs_free(w->scal);
    }
    free_residuals(w->r_orig);
    if (w->stgs && w->stgs->normalize) {
      SCS(free_sol)(w->xys_normalized);
      free_residuals(w->r_normalized);
    }
    if (w->stgs) {
      if (w->stgs->log_csv_filename)
        scs_free((char *)w->stgs->log_csv_filename);
      if (w->stgs->write_data_filename)
        scs_free((char *)w->stgs->write_data_filename);
      scs_free(w->stgs);
    }
    if (w->k) { /* deep copy */
      SCS(free_cone)(w->k);
    }
    if (w->d) { /* deep copy */
      SCS(free_data)(w->d);
    }
    scs_free(w);
  }
}

/* this just calls scs_init, scs_solve, and scs_finish */
scs_int scs(const ScsData *d, const ScsCone *k, const ScsSettings *stgs,
            ScsSolution *sol, ScsInfo *info) {
  scs_int status;
  ScsWork *w = scs_init(d, k, stgs);
  if (w) {
    scs_solve(w, sol, info, stgs->warm_start);
    status = info->status_val;
  } else {
    status = failure(SCS_NULL, d ? d->m : -1, d ? d->n : -1, sol, info,
                     SCS_FAILED, "could not initialize work", "failure");
  }
  scs_finish(w);
  return status;
}
