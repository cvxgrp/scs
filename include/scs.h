#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>

#include "aa.h"
#include "glbopts.h"

/* stores the necessary private workspace, only the linear system solver
 * interacts with this struct */
typedef struct SCS_LIN_SYS_WORK ScsLinSysWork;
typedef struct SCS_ACCEL_WORK ScsAccelWork;
typedef struct SCS_CONE_WORK ScsConeWork;
typedef struct SCS_MATRIX ScsMatrix;

/** struct containing all settings */
typedef struct {
  /* these *cannot* change for multiple runs with the same call to SCS(init) */
  /** whether to heuristically rescale the data before solve */
  scs_int normalize;
  /** initial dual scaling factor (may be updated if adaptive_scale is on) */
  scs_float scale;
  /** whether to adaptively update `scale` */
  scs_int adaptive_scale;
  /** primal constraint scaling factor */
  scs_float rho_x;
  /** maximum iterations to take */
  scs_int max_iters;
  /** absolute convergence tolerance */
  scs_float eps_abs;
  /** relative convergence tolerance */
  scs_float eps_rel;
  /** infeasible convergence tolerance */
  scs_float eps_infeas;
  /** Douglas-Rachford relaxation parameter */
  scs_float alpha;
  /** time limit in secs (can be fractional) */
  scs_float time_limit_secs;
  /** whether to log progress to stdout */
  scs_int verbose;
  /** whether to use warm start (put initial guess in ScsSolution struct) */
  scs_int warm_start;
  /** memory for acceleration */
  scs_int acceleration_lookback;
  /** interval to apply acceleration */
  scs_int acceleration_interval;
  /** string, if set will dump raw prob data to this file */
  const char *write_data_filename;
  /** string, if set will log solve data to this csv file (slows solve a lot) */
  const char *log_csv_filename;
} ScsSettings;

/** struct containing problem data */
typedef struct {
  /** A has `m` rows */
  scs_int m;
  /** A has `n` cols, P has `n` cols and `n` rows */
  scs_int n;
  /** A is supplied in CSC format (size `m` x `n`) */
  ScsMatrix *A;
  /** P is supplied in CSC format, must be upper triangular  (size `n` x `n`) */
  ScsMatrix *P;
  /** dense array for b (size `m`) */
  scs_float *b;
  /** dense array for c (size `n`) */
  scs_float *c;
} ScsData;

/** Cone data. NB: rows of data matrix `A` must be specified in this exact order
 */
typedef struct {
  /** number of linear equality constraints (primal zero, dual free) */
  scs_int z;
  /** number of positive orthant cones */
  scs_int l;
  /** upper box values, `len(bu) = len(bl) = max(bsize-1, 0)` */
  scs_float *bu;
  /** lower box values, `len(bu) = len(bl) = max(bsize-1, 0)` */
  scs_float *bl;
  /** total length of box cone (includes scale `t`) */
  scs_int bsize;
  /** array of second-order cone constraints */
  scs_int *q;
  /** length of SOC array */
  scs_int qsize;
  /** array of semidefinite cone constraints */
  scs_int *s;
  /** length of semidefinite constraints array */
  scs_int ssize;
  /** number of primal exponential cone triples */
  scs_int ep;
  /** number of dual exponential cone triples */
  scs_int ed;
  /** array of power cone params, must be in `[-1, 1]`, negative values are
   * interpreted as specifying the dual cone */
  scs_float *p;
  /** number of (primal and dual) power cone triples */
  scs_int psize;
} ScsCone;

/** contains primal-dual solution arrays */
typedef struct {
  /** primal variable */
  scs_float *x;
  /** dual variable */
  scs_float *y;
  /** slack variable */
  scs_float *s;
} ScsSolution;

/** contains terminating information */
typedef struct {
  /** number of iterations taken */
  scs_int iter;
  /** status string, e.g. 'solved' */
  char status[64];
  /** status as scs_int, defined in glbopts.h */
  scs_int status_val;
  /** number of updates to scale */
  scs_int scale_updates;
  /** primal objective */
  scs_float pobj;
  /** dual objective */
  scs_float dobj;
  /** primal equality residual */
  scs_float res_pri;
  /** dual equality residual */
  scs_float res_dual;
  /** duality gap */
  scs_float gap;
  /** infeasibility cert residual */
  scs_float res_infeas;
  /** unbounded cert residual */
  scs_float res_unbdd_a;
  /** unbounded cert residual */
  scs_float res_unbdd_p;
  /** time taken for setup phase (milliseconds) */
  scs_float setup_time;
  /** time taken for solve phase (milliseconds) */
  scs_float solve_time;
  /** (final) scale parameter */
  scs_float scale;
} ScsInfo;

/* the following structs are not exposed to user */

/* contains normalization variables */
typedef struct {
  scs_float *D, *E; /* for normalization */
  scs_float primal_scale, dual_scale;
} ScsScaling;

/* to hold residual information, *all are unnormalized* */
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

/* workspace for SCS */
typedef struct {
  /* x_prev = x from previous iteration */
  scs_int time_limit_reached; /* set if the time-limit is reached */
  scs_float *u, *v, *u_t, *v_prev, *rsk;
  scs_float *h;                  /* h = [c; b] */
  scs_float *g;                  /* g = (I + M)^{-1} h */
  scs_float *lin_sys_warm_start; /* linear system warm-start (indirect only) */
  scs_float *rho_y_vec; /* vector of rho y parameters (affects cone project) */
  AaWork *accel;        /* struct for acceleration workspace */
  scs_float *b_orig, *c_orig;             /* original b and c vectors */
  scs_float *b_normalized, *c_normalized; /* normalized b and c vectors */
  scs_int m, n;                           /* A has m rows, n cols */
  ScsMatrix *A;                           /* (possibly normalized) A matrix */
  ScsMatrix *P;                           /* (possibly normalized) P matrix */
  ScsLinSysWork *p;            /* struct populated by linear system solver */
  ScsScaling *scal;            /* contains the re-scaling data */
  ScsConeWork *cone_work;      /* workspace for the cone projection step */
  scs_int *cone_boundaries;    /* array with boundaries of cones */
  scs_int cone_boundaries_len; /* total length of cones */
  /* normalized and unnormalized residuals */
  ScsResiduals *r_orig, *r_normalized;
  /* track x,y,s as alg progresses, tau *not* divided out */
  ScsSolution *xys_orig, *xys_normalized;
  /* updating scale params workspace */
  scs_float sum_log_scale_factor;
  scs_int last_scale_update_iter, n_log_scale_factor, scale_updates;
  /* aa norm stat */
  scs_float aa_norm;
  scs_float setup_time; /* time taken for setup phase (milliseconds) */
  scs_float scale;      /* current scale parameter */
  const ScsData *d;
  const ScsCone *k;
  const ScsSettings *stgs; /* contains solver settings specified by user */
} ScsWork;

/*
 * main library API
 */

/**
 * Initialize SCS and allocate memory.
 *
 * All the inputs must be already allocated in memory before calling.
 *
 * It performs:
 * - data and settings validation
 * - problem data scaling
 * - automatic parameters tuning (if enabled)
 * - setup linear system solver:
 *      - direct solver: KKT matrix factorization is performed here
 *      - indirect solver: KKT matrix preconditioning is performed here
 * - solve the first linear system
 *
 *
 * @param  d 		 Problem data
 * @param  k 		 Cone data
 * @param  stgs  SCS solver settings
 * @return       Solver work struct
 */
ScsWork *SCS(init)(const ScsData *d, const ScsCone *k, const ScsSettings *stgs);

/**
 * Solve quadratic cone program initialized by SCS(init).
 *
 * @param  w     Workspace allocated by init
 * @param  sol 	 Solver solution struct, will contain solution at termination
 * @param  info  Solver info reporting
 * @return       Flag containing problem status (see \a glbopts.h)
 */
scs_int SCS(solve)(ScsWork *w, ScsSolution *sol, ScsInfo *info);

/**
 * Clean up allocated SCS workspace.
 *
 * @param  w  Workspace allocated by init, will be deallocated.
 */
void SCS(finish)(ScsWork *w);

/**
 * Solve quadratic cone program defined by data in d and cone k.
 *
 * All the inputs must already be allocated in memory before calling.
 *
 * @param  d 		 Problem data
 * @param  k 		 Cone data
 * @param  stgs  SCS solver settings
 * @param  sol   Solution will be stored here
 * @param  info  Information about the solve will be stored here
 * @return       Flag that determines solve type (see \a glbopts.h)
 */
scs_int scs(const ScsData *d, const ScsCone *k, const ScsSettings *stgs,
            ScsSolution *sol, ScsInfo *info);

/**
 * Helper function to set all settings to default values (see \a glbopts.h)
 *
 * @param  stgs  Settings struct that will be populated
 */
void SCS(set_default_settings)(ScsSettings *stgs);

const char *SCS(version)(void);
size_t SCS(sizeof_int)(void);
size_t SCS(sizeof_float)(void);

#ifdef __cplusplus
}
#endif
#endif
