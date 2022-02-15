/*
 * Define ScsWork and related internal-only structs (not part of external API).
 */

#ifndef SCS_WORK_H_GUARD
#define SCS_WORK_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

/** Contains normalization variables. */
typedef struct {
  scs_float *D, *E; /* for normalization */
  scs_int m;        /* Length of D */
  scs_int n;        /* Length of E */
  scs_float primal_scale, dual_scale;
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
  scs_float *b_orig, *c_orig; /* original unnormalized b and c vectors */
  AaWork *accel;              /* struct for acceleration workspace */
  ScsData *d;                 /* Problem data deep copy NORMALIZED */
  ScsCone *k;                 /* Problem cone deep copy */
  ScsSettings *stgs;          /* contains solver settings specified by user */
  ScsLinSysWork *p;           /* struct populated by linear system solver */
  ScsScaling *scal;           /* contains the re-scaling data */
  ScsConeWork *cone_work;     /* workspace for the cone projection step */
  /* normalized and unnormalized residuals */
  ScsResiduals *r_orig, *r_normalized;
  /* track x,y,s as alg progresses, tau *not* divided out */
  ScsSolution *xys_orig, *xys_normalized;
  /* Scale updating workspace */
  scs_float sum_log_scale_factor;
  scs_int last_scale_update_iter, n_log_scale_factor, scale_updates;
  /* AA stats */
  scs_float aa_norm;
  scs_int rejected_accel_steps, accepted_accel_steps;
};

#ifdef __cplusplus
}
#endif
#endif
