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
  /* scratch, n + m each, nothing kept in them between calls: kkt_res and
   * kkt_corr, the right-hand side [c; -b] and the refined candidate of the
   * solve for g (update_work_cache); kkt_corr also the plain solve that
   * decides at a check of the per-iteration refinement while it is on
   * (solve_direct, solve_indirect), and at a decision the right-hand side
   * whose terms refine_next_on sizes, with kkt_res as its scratch */
  scs_float *kkt_res, *kkt_corr;
  /* per-iteration refinement, reset at the start of each solve (see
   * "Accuracy of each iteration's solve" in glbopts.h): whether it is on,
   * how many residual checks this solve has made, how many iterations were
   * solved with refinement on (direct solvers, solve_direct: the previous
   * solution passed as the warm start, and the linear solver solves plainly
   * instead when it is no better a guess than zero; indirect solvers,
   * solve_indirect: with refinement steps allowed, taken only where the CG
   * solve missed its tolerance), how often the safeguard replaced the
   * warm-started solution by the plain solve (direct solvers), all read by
   * the tests, and the tol it compares with (set_refine_tol, at each
   * convergence check; -1 before the first) */
  scs_int refine_on, refine_checks, refine_on_iters, refine_safeguards;
  scs_float refine_tol;
  /* direct solvers: the previous iteration's solution (n + m) whenever
   * refinement is on, passed as the warm start of the next iteration's
   * solve. Unused (NULL) in indirect builds. */
  scs_float *refine_prev;
  /* indirect solvers: the warm start of each iteration's CG solve, n + m (CG
   * reads the first n; the rest stays zero). Unused (NULL) in direct
   * builds. */
  scs_float *lin_sys_warm_start;
  /* indirect solvers: the size of SCS's residuals that each iteration's CG
   * tolerance is set from (project_lin_sys), saved where SCS computes the
   * residuals for the algorithm (set_cg_tol_res, at each convergence check
   * and at the end of each solve); 0 before the first solve. Unused in direct
   * builds. */
  scs_float cg_tol_res;
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
  /* scale and adaptive_diag_scale as given to scs_init. The solve overwrites
   * stgs->scale with the adapted scale, and setup sets
   * stgs->adaptive_diag_scale to 0 when adaptive_scale is off, so
   * scs_update_settings compares a caller's settings with these instead */
  scs_float scale_init;
  scs_int adaptive_diag_scale_init;
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
  /* the multipliers the next update would install, size m (scratch) */
  scs_float *scale_mults_next;
  /* floor on the row-profile denominators, see row_rel_res */
  scs_float den_floor;
  /* AA stats */
  scs_float aa_norm;
  scs_int rejected_accel_steps, accepted_accel_steps;
};

#ifdef __cplusplus
}
#endif
#endif
