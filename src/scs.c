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
/* MKL's interface-layer value is a bitmask: bit 0 selects the integer
 * width (LP64/ILP64) and bit 1 is the GNU Fortran calling-convention
 * flag, so MKL_Set_Interface_Layer can legitimately return 2 (LP64|GNU)
 * or 3 (ILP64|GNU) -- e.g. when an MKL-backed NumPy initialized MKL
 * first. Only the width bit affects integer-size safety. */
#define MKL_INTERFACE_GNU 2
/* bit 0 of the layer value selects the integer width (LP64/ILP64) */
#define MKL_INTERFACE_WIDTH_MASK MKL_INTERFACE_ILP64
/* Only libmkl_rt provides this; a static MKL link has nothing to negotiate,
 * so the reference is weak and the check is skipped when the symbol is
 * absent (MSVC links mkl_rt dynamically and keeps the strong reference). */
#if defined(__GNUC__)
int MKL_Set_Interface_Layer(int) __attribute__((weak));
#else
int MKL_Set_Interface_Layer(int);
#endif

static const char *scs_mkl_interface_name(int layer) {
  switch (layer) {
  case MKL_INTERFACE_LP64:
    return "LP64";
  case MKL_INTERFACE_ILP64:
    return "ILP64";
  case MKL_INTERFACE_LP64 | MKL_INTERFACE_GNU:
    return "LP64 (GNU)";
  case MKL_INTERFACE_ILP64 | MKL_INTERFACE_GNU:
    return "ILP64 (GNU)";
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
    int actual;
#if defined(__GNUC__)
    if (!MKL_Set_Interface_Layer) {
      /* static MKL: no runtime interface layer exists to query */
      return 0;
    }
#endif
    actual = MKL_Set_Interface_Layer(expected);
    /* MKL returns -1 for an invalid request; without this guard the width
     * comparison below would misread it ((-1 & 1) == 1, i.e. ILP64) */
    if (actual < 0) {
      scs_printf("MKL_Set_Interface_Layer(%d) failed with error %d.\n",
                 expected, actual);
      return -1;
    }
    /* compare only the integer-width bit: the GNU flag is a Fortran
     * calling-convention variant with the same integer sizes */
    if ((actual & MKL_INTERFACE_WIDTH_MASK) !=
        (expected & MKL_INTERFACE_WIDTH_MASK)) {
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
static void set_info_aa_stats(ScsInfo *info, const AaWork *accel);
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
#if INDIRECT > 0
static scs_int refine_iter(ScsWork *w, scs_int iter);
#else
static scs_int solve_direct(ScsWork *w, scs_int iter);
#endif
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
  if (!r->ax || !r->ax_s || !r->ax_s_btau || !r->px || !r->aty ||
      !r->px_aty_ctau) {
    free_residuals(r);
    return SCS_NULL;
  }
  return r;
}

static void set_info_aa_stats(ScsInfo *info, const AaWork *accel) {
  if (!info) {
    return;
  }
  memset(&info->aa_stats, 0, sizeof(info->aa_stats));
  info->aa_stats.last_aa_norm = NAN;
  if (accel) {
    info->aa_stats = aa_get_stats(accel);
  }
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
      if (sol->x) {
        SCS(scale_array)(sol->x, NAN, n);
      }
    }
    if (m > 0) {
      if (!sol->y) {
        sol->y = (scs_float *)scs_calloc(m, sizeof(scs_float));
      }
      if (sol->y) {
        SCS(scale_array)(sol->y, NAN, m);
      }
      if (!sol->s) {
        sol->s = (scs_float *)scs_calloc(m, sizeof(scs_float));
      }
      if (sol->s) {
        SCS(scale_array)(sol->s, NAN, m);
      }
    }
  }
}

static scs_int failure(ScsWork *w, scs_int m, scs_int n, ScsSolution *sol,
                       ScsInfo *info, scs_int stint, const char *msg,
                       const char *ststr) {
  scs_int status = stint;
  SCS(close_csv_log_file)(w);
  populate_on_failure(m, n, sol, info, status, ststr);
  set_info_aa_stats(info, w ? w->accel : SCS_NULL);
  scs_printf("Failure:%s\n", msg);
  scs_end_interrupt_listener();
  return status;
}

/* ========================= Validation ============================== */

#if NO_VALIDATE == 0
/* Checks the values of the settings alone; scs_init and scs_update_settings
 * both use it, so a value is accepted by one exactly when by the other. */
static scs_int validate_settings(const ScsSettings *stgs);

static scs_int validate(const ScsData *d, const ScsCone *k,
                        const ScsSettings *stgs) {
  if (d->m <= 0 || d->n <= 0) {
    scs_printf("m and n must both be greater than 0; m = %li, n = %li\n",
               (long)d->m, (long)d->n);
    return -1;
  }
  if (d->A && (d->A->m != d->m || d->A->n != d->n)) {
    scs_printf("A dimensions = (%li, %li), inconsistent with m = %li, "
               "n = %li\n",
               (long)d->A->m, (long)d->A->n, (long)d->m, (long)d->n);
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
  return validate_settings(stgs);
}

static scs_int validate_settings(const ScsSettings *stgs) {
  if (stgs->max_iters <= 0) {
    scs_printf("max_iters must be positive\n");
    return -1;
  }
  if (!isfinite(stgs->eps_abs) || stgs->eps_abs < 0) {
    scs_printf("eps_abs tolerance must be a nonnegative finite number\n");
    return -1;
  }
  if (!isfinite(stgs->eps_rel) || stgs->eps_rel < 0) {
    scs_printf("eps_rel tolerance must be a nonnegative finite number\n");
    return -1;
  }
  if (!isfinite(stgs->eps_infeas) || stgs->eps_infeas < 0) {
    scs_printf("eps_infeas tolerance must be a nonnegative finite number\n");
    return -1;
  }
  if (!isfinite(stgs->alpha) || stgs->alpha <= 0 || stgs->alpha >= 2) {
    scs_printf("alpha must be in (0,2)\n");
    return -1;
  }
  if (!isfinite(stgs->rho_x) || stgs->rho_x <= 0) {
    scs_printf("rho_x must be a positive finite number (1e-3 works well).\n");
    return -1;
  }
  if (stgs->adaptive_diag_scale < 0 || stgs->adaptive_diag_scale > 1) {
    scs_printf("adaptive_diag_scale must be 0 (off) or 1 (on).\n");
    return -1;
  }
  if (!isfinite(stgs->scale) || stgs->scale <= 0) {
    scs_printf("scale must be a positive finite number (1 works well).\n");
    return -1;
  }
  if (!isfinite(stgs->time_limit_secs) || stgs->time_limit_secs < 0) {
    scs_printf("time_limit_secs must be a nonnegative finite number.\n");
    return -1;
  }
  if (stgs->acceleration_interval <= 0) {
    scs_printf("acceleration_interval must be positive (10 works well).\n");
    return -1;
  }
  if (stgs->acceleration_lookback < 0) {
    scs_printf("acceleration_lookback must be nonnegative "
               "(use acceleration_type_1=0 for type-II AA).\n");
    return -1;
  }
  if (!isfinite(stgs->acceleration_regularization)) {
    /* Sign-encoded modes per include/aa.h: positive = scaled by
     * ||A||_F ||Y||_F, negative = pinned absolute value, zero = off. */
    scs_printf("acceleration_regularization must be a finite number.\n");
    return -1;
  }
  if (!isfinite(stgs->acceleration_relaxation) ||
      stgs->acceleration_relaxation < 0 ||
      stgs->acceleration_relaxation > 2) {
    scs_printf("acceleration_relaxation must be in [0, 2].\n");
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
  /* Certificate residuals divide the violation by the certificate's
   * objective value (-c'x or -b'y), which bakes the data magnitudes into
   * the effective tolerance: with e.g. ||c|| ~ 1e6 the unboundedness test
   * is six orders looser than intended and produces false certificates
   * (netlib agg/grow families). Tighten by the data-norm ratio so the
   * test is invariant to rescaling of c (resp. b) vs A; only ever
   * tighten. Additionally require the certificate to hold on consecutive
   * checks to reject transient (e.g. accelerated) iterates. */
  {
    scs_float tighten_unbdd =
        MAX(1., SAFEDIV_POS(w->nm_c_orig, w->nm_a_orig));
    scs_float tighten_infeas =
        MAX(1., SAFEDIV_POS(w->nm_b_orig, w->nm_a_orig));
    if (isless(r->res_unbdd_a * tighten_unbdd, eps_infeas) &&
        isless(r->res_unbdd_p, eps_infeas)) {
      if (++w->unbdd_cert_streak >= CERT_PERSISTENCE_CHECKS) {
        return SCS_UNBOUNDED;
      }
    } else {
      w->unbdd_cert_streak = 0;
    }
    if (isless(r->res_infeas * tighten_infeas, eps_infeas)) {
      if (++w->infeas_cert_streak >= CERT_PERSISTENCE_CHECKS) {
        return SCS_INFEASIBLE;
      }
    } else {
      w->infeas_cert_streak = 0;
    }
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

static scs_float root_plus_from_coeffs(scs_float a, scs_float b, scs_float c) {
  scs_float rad, sqrt_rad, q;
  if (!isfinite(a) || !isfinite(b) || !isfinite(c) || a <= 0.) {
    return NAN;
  }
  rad = b * b - 4 * a * c;
  if (!isfinite(rad)) {
    return NAN;
  }
  if (rad < 0.) {
    /* Preserve the historical repeated-root fallback for no-real-root cases. */
    return -b / (2 * a);
  }
  sqrt_rad = SQRTF(rad);
  if (b <= 0.) {
    return (-b + sqrt_rad) / (2 * a);
  }
  q = -0.5 * (b + sqrt_rad);
  return q != 0. ? c / q : 0.;
}

static scs_float root_plus(ScsWork *w, scs_float *p, scs_float *mu,
                           scs_float eta) {
  /* Compute all five weighted dot products (g'Rg, mu'Rg, p'Rg, p'Rp, p'Rmu)
   * in a single pass over diag_r to minimise memory traffic. */
  scs_int i, nm = w->d->n + w->d->m;
  scs_float gg = 0., mug = 0., pg = 0., pp = 0., pmu = 0.;
  scs_float a, b, c, tau_scale = w->diag_r[nm];
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
  return root_plus_from_coeffs(a, b, c);
}

/* b = [R_x v_x; -R_y v_y], the right-hand side of this iteration's solve
 * (n + m entries) */
static void iter_rhs(const ScsWork *w, scs_float *b) {
  scs_int i, n = w->d->n, nm = n + w->d->m;
  /* Copy and scale in one pass, eliminating the intermediate memcpy. */
  for (i = 0; i < n; ++i) {
    b[i] = w->v[i] * w->diag_r[i];
  }
  for (i = n; i < nm; ++i) {
    b[i] = -w->v[i] * w->diag_r[i];
  }
}

/* status != 0 indicates failure */
static scs_int project_lin_sys(ScsWork *w, scs_int iter) {
  scs_int n = w->d->n, m = w->d->m, l = n + m + 1, status;
#if INDIRECT > 0
  scs_float *warm_start = SCS_NULL;
  scs_float tol = -1.0;
  scs_float nm_ax_s_btau, nm_px_aty_ctau, nm_ws;
  iter_rhs(w, w->u_t);
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
  status = scs_solve_lin_sys(w->p, w->u_t, warm_start, tol);
  if (status == 0) {
    status = refine_iter(w, iter);
  }
#else
  status = solve_direct(w, iter);
#endif
  w->u_t[l - 1] = w->v[l - 1];
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
  set_info_aa_stats(info, w->accel);
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
  if (w->scale_mults) {
    /* per-row multipliers act like a per-row scale: R_y_i = base_i / f_i.
     * Multipliers are kept in the form their cone block needs (enforced
     * at update time) so no SCS(enforce_cone_boundaries) needed here. */
    for (i = 0; i < w->d->m; ++i) {
      w->diag_r[w->d->n + i] /= w->scale_mults[i];
    }
  }
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
  /* before either is changed below or by the solve, see scs_work.h */
  w->scale_init = w->stgs->scale;
  w->adaptive_diag_scale_init = w->stgs->adaptive_diag_scale;

  /* record the original (pre-normalization) magnitude of A for the
   * data-scale-invariant certificate checks */
  w->nm_a_orig = SCS(norm_inf)(w->d->A->x, w->d->A->p[w->d->A->n]);

  /* allocate workspace: */
  w->u = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->u_t = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->v_prev = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->rsk = (scs_float *)scs_calloc(l, sizeof(scs_float));
  w->h = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->g = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->kkt_res = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
  w->kkt_corr = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
#if INDIRECT == 0
  w->refine_prev = (scs_float *)scs_calloc((l - 1), sizeof(scs_float));
#endif
  w->lin_sys_warm_start = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
  w->diag_r = (scs_float *)scs_calloc(l, sizeof(scs_float));
  /* x,y,s struct */
  w->xys_orig = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  if (!w->xys_orig) {
    scs_printf("ERROR: work memory allocation failure\n");
    scs_finish(w);
    return SCS_NULL;
  }
  w->xys_orig->x = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
  w->xys_orig->s = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  w->xys_orig->y = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  w->r_orig = init_residuals(w->d);
  w->b_orig = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
  w->c_orig = (scs_float *)scs_calloc(w->d->n, sizeof(scs_float));
  if (w->stgs->adaptive_diag_scale) {
    if (!w->stgs->adaptive_scale) {
      /* silently disable: diag scaling rides the adaptive-scale update
       * machinery, and since it defaults on, users who only turn off
       * adaptive_scale should not see a warning */
      w->stgs->adaptive_diag_scale = 0;
    } else {
      scs_int j;
      w->scale_mults = (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
      w->scale_mults_next =
          (scs_float *)scs_calloc(w->d->m, sizeof(scs_float));
      if (!w->scale_mults || !w->scale_mults_next) {
        scs_printf("ERROR: work memory allocation failure\n");
        scs_finish(w);
        return SCS_NULL;
      }
      for (j = 0; j < w->d->m; ++j) {
        w->scale_mults[j] = 1.0;
      }
    }
  }

  if (!w->u || !w->u_t || !w->v || !w->v_prev || !w->rsk || !w->h || !w->g ||
      !w->kkt_res || !w->kkt_corr ||
#if INDIRECT == 0
      !w->refine_prev ||
#endif
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
  if (w->stgs->verbose && w->cone_work->psd_nthreads > 1) {
    scs_printf("\t  psd cones projected in parallel by %li threads\n",
               (long)w->cone_work->psd_nthreads);
  }
  set_diag_r(w);

  if (w->stgs->normalize) {
    w->xys_normalized = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
    if (!w->xys_normalized) {
      scs_printf("ERROR: normalized work memory allocation failure\n");
      scs_finish(w);
      return SCS_NULL;
    }
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
    w->scal = SCS(normalize_a_p)(w->d->P, w->d->A, w->d->b, w->d->c,
                                 w->cone_work);
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
#if INDIRECT == 0
  /* does the linear solver compute the refinement's products with its own
   * K (linsys.h)? kkt_res and kkt_corr are zero here, so this is 0 - K 0.
   * "Not provided": SCS computes them. Any other negative value is an error
   * of the linear solver, and setup fails, as when scs_init_lin_sys_work
   * fails: SCS does not switch to its own products when the solver's fail,
   * so a solve does not depend on a transient error. A NaN (K holds an inf)
   * counts as provided, and the solve deals with it as with any NaN. */
  {
    scs_float probe = scs_lin_sys_residual(w->p, w->kkt_res, w->kkt_corr);
    if (probe == SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED) {
      w->refine_hook = 0;
    } else if (probe < 0) {
      scs_printf("ERROR: the linear solver's products with the KKT matrix "
                 "failed at setup\n");
      scs_finish(w);
      return SCS_NULL;
    } else {
      w->refine_hook = 1;
    }
  }
#endif
  if (w->stgs->acceleration_lookback) {
    /* min_len = mem preserves the previous FILL_MEMORY_BEFORE_SOLVE
       gate — AA holds off until the sliding window is full. */
    if (!(w->accel = aa_init(l, w->stgs->acceleration_lookback,
                             w->stgs->acceleration_lookback,
                             w->stgs->acceleration_type_1,
                             w->stgs->acceleration_regularization,
                             w->stgs->acceleration_relaxation,
                             AA_SAFEGUARD_FACTOR,
                             AA_MAX_WEIGHT_NORM, AA_IR_MAX_STEPS,
                             VERBOSITY))) {
      if (w->stgs->verbose) {
        scs_printf("WARN: aa_init returned NULL, no acceleration applied.\n");
      }
    }
  } else {
    w->accel = SCS_NULL;
  }
  return w;
}

/* r = [c; -b] - K g with K = [R_x + P, A'; A, -R_y], from the (normalized)
 * problem data, so it is the same for every linear solver. Returns
 * ||r||_inf divided by the largest inf-norm among [c; -b], R g, A' g_y,
 * P g_x + A' g_y and A g_x: the size of the terms whose rounding errors the
 * residual is compared with (P g_x on its own is at most the sum of the two
 * that contain it). */
static scs_float g_residual(const ScsWork *w, const scs_float *g,
                            scs_float *r) {
  scs_int i, n = w->d->n, m = w->d->m;
  scs_float *rx = r, *ry = &(r[n]), rg, scale;
  const scs_float *gx = g, *gy = &(g[n]);
  memset(r, 0, (n + m) * sizeof(scs_float));
  SCS(accum_by_atrans)(w->d->A, gy, rx);
  scale = MAX(SCS(norm_inf)(w->d->c, n), SCS(norm_inf)(w->d->b, m));
  scale = MAX(scale, SCS(norm_inf)(rx, n));
  if (w->d->P) {
    SCS(accum_by_p)(w->d->P, gx, rx);
    scale = MAX(scale, SCS(norm_inf)(rx, n));
  }
  SCS(accum_by_a)(w->d->A, gx, ry);
  scale = MAX(scale, SCS(norm_inf)(ry, m));
  for (i = 0; i < n; ++i) {
    rg = w->diag_r[i] * gx[i];
    scale = MAX(scale, ABS(rg));
    rx[i] = w->d->c[i] - (rx[i] + rg);
  }
  for (i = 0; i < m; ++i) {
    rg = w->diag_r[n + i] * gy[i];
    scale = MAX(scale, ABS(rg));
    ry[i] = -w->d->b[i] - (ry[i] - rg);
  }
  return SCS(norm_inf)(r, n + m) / MAX(scale, _DIV_EPS_TOL);
}

/* Iterative refinement of g, see G_SOLVE_TOL in glbopts.h. */
static void refine_g(ScsWork *w) {
  scs_int k, l = w->d->n + w->d->m;
  scs_float res, res_new;
  res = g_residual(w, w->g, w->kkt_res);
  for (k = 0; k < G_REFINE_MAX_STEPS && res > G_SOLVE_TOL; ++k) {
    /* the candidate g + d is formed in kkt_corr and kept only if it is
     * better, so a rejected step leaves g exactly as it was */
    memcpy(w->kkt_corr, w->kkt_res, l * sizeof(scs_float));
    scs_solve_lin_sys(w->p, w->kkt_corr, SCS_NULL, CG_BEST_TOL);
    SCS(add_scaled_array)(w->kkt_corr, w->g, l, 1.);
    res_new = g_residual(w, w->kkt_corr, w->kkt_res);
    if (!(res_new <= res)) {
      break; /* the step made it worse (or is not finite): keep g */
    }
    memcpy(w->g, w->kkt_corr, l * sizeof(scs_float));
    if (res_new > 0.5 * res) {
      break; /* no longer contracting */
    }
    res = res_new;
  }
}

/* r = rhs - K u with K = [R_x + P, A'; A, -R_y] and rhs = [R_x v_x; -R_y v_y],
 * the right-hand side of this iteration's solve in project_lin_sys, from the
 * (normalized) problem data, so it is the same for every linear solver.
 * Returns ||r||_inf, and ||rhs||_inf in *nrhs unless nrhs is SCS_NULL. rhs
 * is rebuilt from v and diag_r rather than copied (v does not change in
 * project_lin_sys), in its own statement so that it is rounded as it is in
 * iter_rhs. r may not alias u. */
static scs_float iter_residual(const ScsWork *w, const scs_float *u,
                               scs_float *r, scs_float *nrhs) {
  scs_int i, n = w->d->n, m = w->d->m;
  const scs_float *v = w->v, *dr = w->diag_r;
  scs_float rhs, nb = 0.;
  memset(r, 0, (n + m) * sizeof(scs_float));
  SCS(accum_by_atrans)(w->d->A, &(u[n]), r);
  if (w->d->P) {
    SCS(accum_by_p)(w->d->P, u, r);
  }
  SCS(accum_by_a)(w->d->A, u, &(r[n]));
  for (i = 0; i < n; ++i) {
    rhs = v[i] * dr[i];
    nb = MAX(nb, ABS(rhs));
    r[i] = rhs - (r[i] + dr[i] * u[i]);
  }
  for (i = n; i < n + m; ++i) {
    rhs = -v[i] * dr[i];
    nb = MAX(nb, ABS(rhs));
    r[i] = rhs - (r[i] - dr[i] * u[i]);
  }
  if (nrhs) {
    *nrhs = nb;
  }
  return SCS(norm_inf)(r, n + m);
}

/* Sets the tol that the per-iteration refinement (refine_iter for indirect
 * solvers, solve_direct for direct ones) compares with, at a convergence
 * check, right after the residuals are computed: the size of the residuals
 * SCS is driving to zero (normalized, tau not divided out), but not below
 * REFINE_TOL_FLOOR:
 * one of them can be zero (an infeasibility or unboundedness certificate),
 * and refinement would then chase rounding errors on every remaining
 * iteration. It is kept here rather than read from r_normalized when the
 * refinement runs because logging to CSV (log_csv_filename) recomputes
 * r_normalized every iteration, and logging must not change the iterates. */
static void set_refine_tol(ScsWork *w) {
  scs_float tol = MIN(SCS(norm_inf)(w->r_normalized->ax_s_btau, w->d->m),
                      SCS(norm_inf)(w->r_normalized->px_aty_ctau, w->d->n));
  w->refine_tol = MAX(tol, REFINE_TOL_FLOOR);
}

#if INDIRECT > 0
/* Indirect solvers: refinement of this iteration's solve, see
 * REFINE_CHECK_INTERVAL, REFINE_INDIRECT_ON_FACTOR and
 * REFINE_INDIRECT_OFF_FACTOR in glbopts.h. Called right after the solve in
 * project_lin_sys, while u_t[:n+m] holds its solution. Returns the status of
 * the refinement solve (0 when there is none). */
static scs_int refine_iter(ScsWork *w, scs_int iter) {
  scs_int status, n = w->d->n, m = w->d->m;
  scs_float res;
  if (!w->refine_on && iter % REFINE_CHECK_INTERVAL != 0) {
    return 0;
  }
  if (w->refine_tol < 0) {
    return 0; /* no convergence check yet, so no residuals to compare with */
  }
  res = iter_residual(w, w->u_t, w->kkt_res, SCS_NULL);
  w->refine_checks++;
  /* on above REFINE_INDIRECT_ON_FACTOR * tol, off only at or below
   * REFINE_INDIRECT_OFF_FACTOR * tol: the gap keeps it from switching off
   * and on again */
  w->refine_on = res > (w->refine_on ? REFINE_INDIRECT_OFF_FACTOR
                                     : REFINE_INDIRECT_ON_FACTOR) *
                           w->refine_tol;
  if (!w->refine_on) {
    return 0;
  }
  /* u_t += K^{-1} r: one solve, cold, to the CG solvers' best accuracy */
  status = scs_solve_lin_sys(w->p, w->kkt_res, SCS_NULL, CG_BEST_TOL);
  if (status == 0) {
    SCS(add_scaled_array)(w->u_t, w->kkt_res, n + m, 1.);
    w->refine_steps++;
  }
  return status;
}

#else

/* Direct solvers: u_t[:n+m], the solution of this iteration's solve, becomes
 * the guess of the next iteration's correction form (kept by the linear
 * solver when it computes the products, else in refine_prev). Called after
 * that solve with no other solve in between, as scs_lin_sys_keep requires.
 * Whenever refinement is on after an iteration, this has been called in it,
 * so the guess is always valid while on. */
static void keep_prev(ScsWork *w) {
  if (w->refine_hook) {
    scs_lin_sys_keep(w->p, w->u_t);
  } else {
    memcpy(w->refine_prev, w->u_t, (w->d->n + w->d->m) * sizeof(scs_float));
  }
}

/* Direct solvers: ||rhs - K u||_inf for this iteration's rhs, with the linear
 * solver's product when it has one (rhs built in kkt_res), else SCS's
 * (iter_residual, r in kkt_res), and ||rhs||_inf in *nrhs unless nrhs is
 * SCS_NULL. u may not be kkt_res. A negative value means the linear solver
 * failed. */
static scs_float check_residual(ScsWork *w, const scs_float *u,
                                scs_float *nrhs) {
  if (w->refine_hook) {
    iter_rhs(w, w->kkt_res);
    if (nrhs) {
      *nrhs = SCS(norm_inf)(w->kkt_res, w->d->n + w->d->m);
    }
    return scs_lin_sys_residual(w->p, w->kkt_res, u);
  }
  return iter_residual(w, u, w->kkt_res, nrhs);
}

/* Direct solvers: this iteration's solve in correction form from the
 * previous iteration's solution u0 = refine_prev: u_t = u0 + d with
 * K d = rhs - K u0. Returns 0, leaving u_t untouched, when the guess is no
 * better than zero (||rhs - K u0|| >= ||rhs||, e.g. after a scale update
 * changed K and rhs, or not finite): the caller then solves plainly.
 * Otherwise returns 1 with the solve's status in *status. */
static scs_int correction_solve(ScsWork *w, scs_int *status) {
  scs_int nm = w->d->n + w->d->m;
  scs_float nr, nrhs;
  if (w->refine_hook) {
    /* the linear solver computes r = rhs - K u0 with its own K and guess
     * and keeps r; u_t holds rhs, for the plain solve if r is not used */
    iter_rhs(w, w->u_t);
    nrhs = SCS(norm_inf)(w->u_t, nm);
    nr = scs_lin_sys_residual(w->p, w->u_t, SCS_NULL);
    if (nr < 0) {
      *status = -1; /* the linear solver failed */
      return 1;
    }
    if (!(nr < nrhs)) {
      return 0;
    }
    *status = scs_lin_sys_correct(w->p, w->u_t); /* u0 + K^{-1} r */
  } else {
    /* r = rhs - K u0, computed into u_t, which the solve then overwrites */
    nr = iter_residual(w, w->refine_prev, w->u_t, &nrhs);
    if (!(nr < nrhs)) {
      return 0;
    }
    *status = scs_solve_lin_sys(w->p, w->u_t, SCS_NULL, -1.);
    if (*status == 0) {
      SCS(add_scaled_array)(w->u_t, w->refine_prev, nm, 1.); /* u0 + d */
    }
  }
  if (*status == 0) {
    w->refine_corrections++;
  }
  return 1;
}

/* Direct solvers: this iteration's solve, with refinement in correction
 * form, see "Accuracy of each iteration's solve" in glbopts.h. Leaves the
 * solution in u_t[:n+m]. Every decision below is a function of the iterates
 * and the residuals, never of timings. While off: a plain solve, and at every
 * REFINE_CHECK_INTERVAL-th iteration its residual is compared with tol (on
 * above REFINE_DIRECT_ON_FACTOR * tol). While on: the correction form from
 * the previous iteration's solution (plain when that guess is no better than
 * zero); at a check, a plain solve of the same right-hand side, into
 * scratch, decides whether to stay on (off at or below
 * REFINE_DIRECT_OFF_FACTOR * tol), and the iterate stays the correction-form
 * solution unless its residual is clearly larger than the plain solve's (see
 * REFINE_DIRECT_SAFEGUARD_FACTOR). Returns the status of the last solve (0:
 * success). */
static scs_int solve_direct(ScsWork *w, scs_int iter) {
  scs_int status = 0, corrected = 0, was_on = w->refine_on;
  scs_int nm = w->d->n + w->d->m;
  scs_float res, res_corr, nrhs;
  if (was_on) {
    corrected = correction_solve(w, &status);
  }
  if (!corrected) {
    iter_rhs(w, w->u_t);
    status = scs_solve_lin_sys(w->p, w->u_t, SCS_NULL, -1.);
  }
  if (status != 0) {
    return status;
  }
  if (iter % REFINE_CHECK_INTERVAL != 0 || w->refine_tol < 0) {
    /* no check (none before the first convergence check, which sets tol) */
    if (was_on) {
      keep_prev(w);
    }
    return 0;
  }
  if (corrected) {
    /* the decision is made on a plain solve of this right-hand side, in
     * kkt_corr; the iterate stays the correction-form solution. (Using the
     * plain solve as the iterate instead mixes an inexact solve into the
     * accurate ones every REFINE_CHECK_INTERVAL iterations, which took up to
     * 31 times the iterations on a badly scaled QP.) */
    keep_prev(w);
    iter_rhs(w, w->kkt_corr);
    status = scs_solve_lin_sys(w->p, w->kkt_corr, SCS_NULL, -1.);
    if (status != 0) {
      return status;
    }
    res = check_residual(w, w->kkt_corr, SCS_NULL);
    if (res < 0) {
      return -1; /* the linear solver failed */
    }
    /* Safeguard, see REFINE_DIRECT_SAFEGUARD_FACTOR: the correction form's
     * own residual (one more multiply per check while on); if it is clearly
     * larger than the plain solve's, or NaN, the iterate and the next guess
     * become the plain solve. */
    res_corr = check_residual(w, w->u_t, &nrhs);
    if (res_corr < 0) {
      return -1; /* the linear solver failed */
    }
    if (res_corr != res_corr ||
        (res == res && res_corr > REFINE_DIRECT_SAFEGUARD_FACTOR * res &&
         res_corr > REFINE_DIRECT_SAFEGUARD_FLOOR * nrhs)) {
      memcpy(w->u_t, w->kkt_corr, nm * sizeof(scs_float));
      /* the plain solve was the linear solver's latest solve, as
       * scs_lin_sys_keep requires */
      keep_prev(w);
    }
  } else {
    res = check_residual(w, w->u_t, SCS_NULL);
    if (res < 0) {
      return -1; /* the linear solver failed (NaN iterates give NaN: off) */
    }
  }
  w->refine_checks++;
  /* on above REFINE_DIRECT_ON_FACTOR * tol, off only at or below
   * REFINE_DIRECT_OFF_FACTOR * tol: the gap keeps it from switching off and
   * on again */
  w->refine_on = res > (was_on ? REFINE_DIRECT_OFF_FACTOR
                               : REFINE_DIRECT_ON_FACTOR) *
                           w->refine_tol;
  if (w->refine_on && !corrected) {
    keep_prev(w);
  }
  return 0;
}
#endif

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
  refine_g(w);
}

/* Reset quantities specific to current solve */
static void reset_tracking(ScsWork *w) {
  w->last_scale_update_iter = 0;
  w->sum_log_scale_factor = 0.;
  w->n_log_scale_factor = 0;
  w->scale_updates = 0;
  w->time_limit_reached = 0;
  /* certificate persistence */
  w->infeas_cert_streak = 0;
  w->unbdd_cert_streak = 0;
  /* Acceleration */
  w->rejected_accel_steps = 0;
  w->accepted_accel_steps = 0;
  w->aa_norm = 0.;
  /* Need this to force residual calc if previous solve solved at iter 0 */
  w->r_normalized->last_iter = -1;
  w->r_orig->last_iter = -1;
  /* per-iteration refinement starts off and checks at its own interval, from
   * the first convergence check on */
  w->refine_on = 0;
  w->refine_checks = 0;
  w->refine_steps = 0;
  w->refine_corrections = 0;
  w->refine_tol = -1.;
}

static scs_int update_work(ScsWork *w, ScsSolution *sol) {
  reset_tracking(w);

  if (w->stgs->warm_start) {
    warm_start_vars(w, sol);
  } else {
    cold_start_vars(w);
  }

  /* the cached spectral warm starts describe the iterates of the previous
   * solve, which the (re)start above has just discarded */
  SCS(reset_cone_cache)(w->cone_work);

  update_work_cache(w);
  return 0;
}

/* will update if the factor is outside of range */
static scs_int should_update_r(scs_float factor) {
  return (factor > SQRTF(10.) || factor < 1. / SQRTF(10.));
}

/* Relative primal residual of constraint row i in the normalized space
 * (where diag_r acts). All quantities carry tau consistently. */
/* Relative primal residual of row i, used as the per-row profile.
 *
 * The denominator carries a floor at DEN_FLOOR_FRAC of the block's rms
 * denominator. Without it, a row whose b_i is zero divides only by
 * max(|(Ax)_i|, |s_i|) -- both of which shrink as the iterates converge
 * -- so its ratio inflates even when the row is converging perfectly
 * well, and the profile ends up reporting which terms a row happens to
 * contain rather than how it is converging. Rows with denominators of
 * ordinary size are unaffected; only the collapsing ones are floored. */
static scs_float row_rel_res(const ScsWork *w, scs_int i) {
  const ScsResiduals *r = w->r_normalized;
  scs_float den = MAX(ABS(r->ax[i]), ABS(w->xys_normalized->s[i]));
  den = MAX(den, ABS(w->d->b[i]) * r->tau);
  den = MAX(den, w->den_floor);
  den = MAX(den, _DIV_EPS_TOL);
  return MAX(ABS(r->ax_s_btau[i]), _DIV_EPS_TOL) / den;
}

/* Recompute the row-profile denominator floor from the current iterate. */
static void set_den_floor(ScsWork *w) {
  const ScsResiduals *r = w->r_normalized;
  scs_int i, m = w->d->m;
  scs_float s2 = 0.;
  for (i = 0; i < m; ++i) {
    scs_float d = MAX(ABS(r->ax[i]), ABS(w->xys_normalized->s[i]));
    d = MAX(d, ABS(w->d->b[i]) * r->tau);
    s2 += d * d;
  }
  w->den_floor = DEN_FLOOR_FRAC * SQRTF(s2 / (scs_float)MAX(m, 1));
}

static scs_int update_scale(ScsWork *w, const ScsCone *k, scs_int iter) {
  scs_int i;
  scs_float factor, new_scale, dist_pri, dist_dual;
  scs_float terms_pri, terms_dual, abs_pri, abs_dual, sat_pri, sat_dual, wt;
  const scs_float eps_abs = w->stgs->eps_abs;
  const scs_float eps_rel = w->stgs->eps_rel;

  ScsResiduals *r = w->r_orig;

  scs_float nm_ax = SCALE_NORM(r->ax, w->d->m);
  scs_float nm_s = SCALE_NORM(w->xys_orig->s, w->d->m);
  scs_float nm_px_aty_ctau = SCALE_NORM(r->px_aty_ctau, w->d->n);
  scs_float nm_px = SCALE_NORM(r->px, w->d->n);
  scs_float nm_aty = SCALE_NORM(r->aty, w->d->n);
  scs_float nm_ax_s_btau = SCALE_NORM(r->ax_s_btau, w->d->m);

  scs_int iters_since_last_update = iter - w->last_scale_update_iter;
  if (w->scale_mults) {
    set_den_floor(w);
  }
  /* The scale is moved toward whichever residual is further behind, so each
   * sample is the log of a primal-to-dual ratio. Two ratios are blended.
   *
   * (1) Distance to stopping: each residual over its own threshold in
   * has_converged, with tau multiplied through. At or below one, that side
   * would pass. Balancing these evens out how far each side has to go. */
  /* ||Ax + s - b * tau|| */
  terms_pri = MAX(MAX(nm_ax, nm_s), w->nm_b_orig * r->tau);
  dist_pri = SAFEDIV_POS(nm_ax_s_btau, eps_abs * r->tau + eps_rel * terms_pri);
  /* ||Px + A'y + c * tau|| */
  terms_dual = MAX(MAX(nm_px, nm_aty), w->nm_c_orig * r->tau);
  dist_dual =
      SAFEDIV_POS(nm_px_aty_ctau, eps_abs * r->tau + eps_rel * terms_dual);
  /* clamp to avoid log(0) which would NaN-poison sum_log_scale_factor */
  dist_pri = MAX(dist_pri, _DIV_EPS_TOL);
  dist_dual = MAX(dist_dual, _DIV_EPS_TOL);

  /* (2) The residuals themselves, in the space the iteration runs in
   * (equilibrated when normalize is on), with nothing divided out. */
  abs_pri = MAX(SCALE_NORM(w->r_normalized->ax_s_btau, w->d->m), _DIV_EPS_TOL);
  abs_dual =
      MAX(SCALE_NORM(w->r_normalized->px_aty_ctau, w->d->n), _DIV_EPS_TOL);

  /* A residual as large as the terms it is measured against has no correct
   * digits, and its distance then sits at a ceiling of about 1 / eps_rel
   * whatever the iterate does. With c = 0 and P = 0 the dual is there on
   * every iteration, because the residual A'y is the only term: (1) alone
   * then reads the dual as permanently behind and walks the scale to its
   * floor. sat is 1 at that ceiling and 0 when the residual is negligible
   * against its terms; (1) gets weight wt, the trust in both sides, and (2)
   * gets the rest. */
  sat_pri = MIN(1., eps_rel * dist_pri);
  sat_dual = MIN(1., eps_rel * dist_dual);
  wt = (1. - sat_pri) * (1. - sat_dual);

  /* higher scale makes res_pri go down faster, so increase if res_pri larger */
  w->sum_log_scale_factor += wt * (log(dist_pri) - log(dist_dual)) +
                             (1. - wt) * (log(abs_pri) - log(abs_dual));
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
  scs_int apply_scalar =
      (new_scale != w->stgs->scale) && should_update_r(factor);
  scs_int apply_diag = 0;
  scs_float log_g = 0., g, step, change;
  if (w->stgs->adaptive_diag_scale) {
    /* Residual profiles. The multipliers the update would install are
     * built in full first: a damped step along the row profile, the clamp,
     * and then the form each cone block must take (one value for most
     * non-polyhedral cones, a rank-one fit for PSD blocks). The diagonal
     * adapts when the scalar updates, or when some multiplier would move
     * by more than the update threshold, so that a railed scalar does not
     * freeze the diagonal and a railed multiplier does not keep asking for
     * an update it cannot take.
     *
     * The change is measured on what would be installed, after the cone
     * blocks have been given their form. Measured before that, a single
     * row of a block that wants a large step refactorizes the matrix even
     * though its block then moves by almost nothing: the factorization,
     * the acceleration memory and the cone cache are thrown away as often
     * as every RESCALING_MIN_ITERS iterations for a metric that has not
     * changed. (The fit of a PSD block can also move an entry further than
     * any one row asked for; that change is real and is counted too.)
     * Rows of the zero, nonnegative and box cones are not shared, so for
     * them the two measures are the same. */
    scs_float change_max = 1.0;
    for (i = 0; i < w->d->m; ++i) {
      log_g += log(row_rel_res(w, i));
    }
    log_g /= (scs_float)w->d->m;
    g = exp(log_g);
    for (i = 0; i < w->d->m; ++i) {
      step = POWF(row_rel_res(w, i) / g, DIAG_SCALE_DAMP);
      w->scale_mults_next[i] =
          MIN(MAX(w->scale_mults[i] * step, DIAG_SCALE_MULT_MIN),
              DIAG_SCALE_MULT_MAX);
    }
    SCS(enforce_cone_boundaries)(w->cone_work, w->scale_mults_next,
                                 &SCS(mean), 1);
    for (i = 0; i < w->d->m; ++i) {
      change = w->scale_mults_next[i] / w->scale_mults[i];
      change_max = MAX(change_max, MAX(change, 1. / change));
    }
    apply_diag = apply_scalar || should_update_r(change_max);
  }
  if (!apply_scalar && !apply_diag) {
    return 0;
  }
  {
    scs_int linsys_status;
    w->scale_updates++;
    w->sum_log_scale_factor = 0;
    w->n_log_scale_factor = 0;
    w->last_scale_update_iter = iter;
    if (apply_scalar) {
      w->stgs->scale = new_scale;
    }
    if (apply_diag) {
      /* already clamped and in the form the cone blocks need */
      memcpy(w->scale_mults, w->scale_mults_next,
             w->d->m * sizeof(scs_float));
    }

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
    /* Same reasoning for the cone projections: the spectral cones' inner
     * solvers warm-start from state computed under the metric we just
     * replaced, so drop it rather than let it seed the next iteration's
     * projection (see SCS(reset_cone_cache)). */
    SCS(reset_cone_cache)(w->cone_work);
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

/* Settings that scs_init uses to build the workspace, and why each cannot
 * change once it is built:
 *   normalize        the data in w->d were equilibrated in init_work and the
 *                    normalized iterates allocated (or aliased) there
 *   scale, rho_x     they set diag_r, which is inside the factorized (or
 *                    preconditioned) matrix built by scs_init_lin_sys_work
 *   adaptive_scale   decided in init_work whether the adaptive_diag_scale
 *                    multipliers exist; switching it on or off afterwards
 *                    would leave adaptive_diag_scale meaning something other
 *                    than what was asked for
 *   adaptive_diag_scale  its multipliers are allocated in init_work
 *   acceleration_lookback, _type_1, _regularization, _relaxation
 *                    passed to aa_init, which sized and configured the
 *                    acceleration memory with them
 *   acceleration_interval  read only by the solve loop, but the acceleration
 *                    memory is kept from one solve to the next and holds
 *                    steps taken at the old interval
 *   write_data_filename  the data were written in scs_init; a new name
 *                    would write nothing
 * Changing any of them would need the setup redone, which is what scs_init
 * is for, so the call refuses rather than ignoring the new value or
 * redoing part of the setup behind the caller's back.
 *
 * Exact comparison is intended: the caller passes back the values it gave
 * scs_init, so any difference is a request for a change. */
static scs_int str_differ(const char *a, const char *b) {
  if (!a || !b) {
    return a != b;
  }
  return strcmp(a, b) != 0;
}

static scs_int setup_setting_changed(const char *name, scs_int changed) {
  if (changed) {
    scs_printf("ERROR: scs_update_settings: %s is used to set up the "
               "workspace and cannot be changed after scs_init; it must "
               "equal the value passed to scs_init (create a new workspace "
               "to change it).\n",
               name);
  }
  return changed;
}

static scs_int setup_settings_changed(const ScsWork *w,
                                      const ScsSettings *stgs) {
  const ScsSettings *old = w->stgs;
  scs_int n = 0;
  n += setup_setting_changed("normalize", stgs->normalize != old->normalize);
  n += setup_setting_changed("scale", stgs->scale != w->scale_init);
  n += setup_setting_changed("adaptive_scale",
                             stgs->adaptive_scale != old->adaptive_scale);
  n += setup_setting_changed("adaptive_diag_scale",
                             stgs->adaptive_diag_scale !=
                                 w->adaptive_diag_scale_init);
  n += setup_setting_changed("rho_x", stgs->rho_x != old->rho_x);
  n += setup_setting_changed("acceleration_lookback",
                             stgs->acceleration_lookback !=
                                 old->acceleration_lookback);
  n += setup_setting_changed("acceleration_interval",
                             stgs->acceleration_interval !=
                                 old->acceleration_interval);
  n += setup_setting_changed("acceleration_type_1",
                             stgs->acceleration_type_1 !=
                                 old->acceleration_type_1);
  n += setup_setting_changed("acceleration_regularization",
                             stgs->acceleration_regularization !=
                                 old->acceleration_regularization);
  n += setup_setting_changed("acceleration_relaxation",
                             stgs->acceleration_relaxation !=
                                 old->acceleration_relaxation);
  n += setup_setting_changed("write_data_filename",
                             str_differ(stgs->write_data_filename,
                                        old->write_data_filename));
  return n;
}

scs_int scs_update_settings(ScsWork *w, const ScsSettings *stgs) {
  char *log_csv_filename = SCS_NULL;
  if (!w || !stgs) {
    scs_printf("ERROR: scs_update_settings: missing ScsWork or ScsSettings "
               "input\n");
    return -1;
  }
  if (setup_settings_changed(w, stgs) > 0) {
    scs_printf("ERROR: scs_update_settings: no settings were changed.\n");
    return -1;
  }
#if NO_VALIDATE == 0
  if (validate_settings(stgs) < 0) {
    scs_printf("ERROR: scs_update_settings: no settings were changed.\n");
    return -1;
  }
#endif
  /* copy the one string first: if that fails nothing has been changed */
  if (stgs->log_csv_filename) {
    log_csv_filename = (char *)scs_malloc(strlen(stgs->log_csv_filename) + 1);
    if (!log_csv_filename) {
      scs_printf("ERROR: scs_update_settings: allocation failure, no "
                 "settings were changed.\n");
      return -1;
    }
    strcpy(log_csv_filename, stgs->log_csv_filename);
  }
  /* The settings below are read only by scs_solve (the CSV log is opened at
   * the start of each solve and closed at its end, so no file is open now).
   * The setup settings are left alone: they are equal to the stored ones
   * except scale and adaptive_diag_scale, whose stored values are the ones
   * the solver is using and must stay. */
  w->stgs->max_iters = stgs->max_iters;
  w->stgs->eps_abs = stgs->eps_abs;
  w->stgs->eps_rel = stgs->eps_rel;
  w->stgs->eps_infeas = stgs->eps_infeas;
  w->stgs->alpha = stgs->alpha;
  w->stgs->time_limit_secs = stgs->time_limit_secs;
  w->stgs->verbose = stgs->verbose;
  w->stgs->warm_start = stgs->warm_start;
  if (w->stgs->log_csv_filename) {
    scs_free((char *)w->stgs->log_csv_filename);
  }
  w->stgs->log_csv_filename = log_csv_filename;
  if (w->stgs->verbose) {
    scs_printf("settings updated: eps_abs: %.1e, eps_rel: %.1e, "
               "eps_infeas: %.1e\n\t  alpha: %.2f, max_iters: %i, "
               "time_limit_secs: %.2e\n",
               w->stgs->eps_abs, w->stgs->eps_rel, w->stgs->eps_infeas,
               w->stgs->alpha, (int)w->stgs->max_iters,
               w->stgs->time_limit_secs);
  }
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
  SCS(open_csv_log_file)(w);

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
      set_refine_tol(w);
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
    if (w->log_csv_fout) {
      /* calc residuals every iter if logging to csv */
      populate_residual_struct(w, i);
      SCS(log_data_to_csv)(k, w, i, &solve_timer);
    }
  }

  /* Final logging after full run */
  if (w->log_csv_fout) {
    populate_residual_struct(w, i);
    SCS(log_data_to_csv)(k, w, i, &solve_timer);
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

  SCS(close_csv_log_file)(w);
  scs_end_interrupt_listener();
  return info->status_val;
}

void scs_finish(ScsWork *w) {
  if (w) {
    SCS(close_csv_log_file)(w);
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
    scs_free(w->kkt_res);
    scs_free(w->kkt_corr);
    scs_free(w->refine_prev);
    scs_free(w->b_orig);
    scs_free(w->c_orig);
    scs_free(w->lin_sys_warm_start);
    scs_free(w->diag_r);
    scs_free(w->scale_mults);
    scs_free(w->scale_mults_next);
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
