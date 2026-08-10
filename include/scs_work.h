/*
 * Internal workspace structs used during the solve.
 *
 * ScsWork holds all mutable state for one solve: ADMM iterates, residuals,
 * normalization data, and pointers to the linear system and cone workspaces.
 * ScsScaling holds the diagonal matrices from Ruiz equilibration.
 * ScsResiduals tracks primal/dual residuals and infeasibility certificates.
 *
 * These are internal to SCS -- the public API is in scs.h.
 */

#ifndef SCS_WORK_H_GUARD
#define SCS_WORK_H_GUARD

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include "../linsys/csparse.h"

/** Contains normalization variables. */
typedef struct {
  scs_float *D, *E; /* for normalization */
  scs_int m;        /* Length of D */
  scs_int n;        /* Length of E */
  scs_float primal_scale, dual_scale;
  /* accumulated tau-slot scaling from the stacked-operator equilibration;
   * becomes the sigma applied to b, c in normalize_b_c */
  scs_float tau_scale;
} ScsScaling;

/** Holds residual information. */
typedef struct {
  scs_int last_iter;
  scs_float xt_p_x;     /* x' P x  */
  scs_float xt_p_x_tau; /* x'Px * tau^2 *not* divided out */
  scs_float ctx;
  scs_float ctx_tau; /* tau *not* divided out */
  scs_float bty;
  scs_float bty_tau; /* tau *not* divided out */
  scs_float pobj;    /* primal objective */
  scs_float dobj;    /* dual objective */
  scs_float gap;     /* pobj - dobj */
  scs_float tau;
  scs_float kap;
  scs_float res_pri;
  scs_float res_dual;
  scs_float res_infeas;
  scs_float res_unbdd_p;
  scs_float res_unbdd_a;
  /* tau NOT divided out */
  scs_float *ax, *ax_s, *px, *aty, *ax_s_btau, *px_aty_ctau;
} ScsResiduals;

/** Workspace for SCS. */
struct SCS_WORK {
  /* x_prev = x from previous iteration */
  scs_float setup_time;       /* time taken for setup phase (milliseconds) */
  scs_int time_limit_reached; /* boolean, if the time-limit is reached */
  scs_float *u, *u_t;
  scs_float *v, *v_prev;
  scs_float *rsk;                /* rsk [ r; s; kappa ] */
  scs_float *h;                  /* h = [c; b] */
  scs_float *g;                  /* g = (I + M)^{-1} h */
  scs_float *lin_sys_warm_start; /* linear system warm-start (indirect only) */
  scs_float *diag_r; /* vector of R matrix diagonals (affects cone proj) */
  scs_float *b_orig, *c_orig;     /* original unnormalized b and c vectors */
  scs_float nm_b_orig, nm_c_orig; /* unnormalized NORM(b), NORM(c) */
  scs_float nm_a_orig;            /* unnormalized max abs entry of A */
  /* consecutive checks the infeas/unbdd certificate has passed */
  scs_int infeas_cert_streak, unbdd_cert_streak;
  AaWork *accel;                  /* struct for acceleration workspace */
  ScsData *d;                     /* Problem data deep copy NORMALIZED */
  ScsCone *k;                     /* Problem cone deep copy */
  ScsSettings *stgs;      /* contains solver settings specified by user */
  FILE *log_csv_fout;     /* open CSV log stream for current solve */
  ScsLinSysWork *p;       /* struct populated by linear system solver */
  ScsScaling *scal;       /* contains the re-scaling data */
  ScsConeWork *cone_work; /* workspace for the cone projection step */
  /* normalized and unnormalized residuals */
  ScsResiduals *r_orig, *r_normalized;
  /* track x,y,s as alg progresses, tau *not* divided out */
  ScsSolution *xys_orig, *xys_normalized;
  /* Scale updating workspace */
  scs_float sum_log_scale_factor;
  scs_int last_scale_update_iter, n_log_scale_factor, scale_updates;
  /* per-row scale multipliers on R_y, size m (adaptive_diag_scale >= 1) */
  scs_float *scale_mults;
  /* per-column multipliers on rho_x, size n (adaptive_diag_scale >= 2) */
  scs_float *col_mults;
  /* AA stats */
  scs_float aa_norm;
  scs_int rejected_accel_steps, accepted_accel_steps;
  /* (research) SOC automorphism block metric (stgs->soc_metric):
   * per-block W = r * P(w) with w a unit-determinant boost vector; sw is
   * the cached Jordan square root of w. vals holds the -W KKT entries
   * (upper triangle, column-major per block) registered with the linear
   * system. dir/stat carry the streaming residual-direction tracker. */
  scs_int nsoc, soc_wlen, soc_nnz;
  scs_int *soc_starts, *soc_sizes, *soc_off; /* row offset, size, w offset */
  scs_float *soc_w, *soc_sw, *soc_vals, *soc_dir;
  scs_float *soc_stat; /* per block: ema_top, ema_tot, tau */
  scs_float soc_prev_rel; /* max relative residual at last update event */
  ScsSocMetric soc_sm;
  /* Halpern restart state (active when stgs->restart is set) */
  scs_float *v_anchor;          /* Halpern anchor point, size n+m+1 */
  scs_float *v_avg;             /* running average over current cycle */
  scs_float *restart_scratch;   /* scratch for candidate eval, size 4l */
  scs_int avg_count;            /* iterates accumulated into v_avg */
  scs_int restart_inner;        /* iterations since last restart */
  scs_int last_restart_iter;    /* iteration of last restart */
  scs_float restart_merit_last; /* merit at last restart */
  scs_float restart_merit_prev; /* merit at previous restart check */
  scs_int restarts;             /* total restarts this solve */
};

#ifdef __cplusplus
}
#endif
#endif
