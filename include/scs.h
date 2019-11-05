#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include "glbopts.h"
#include "aa.h"

/* private data structs (that you define) containing any necessary data to solve
 * linear system, etc. this defines the matrix A, only the linear system solver
 * interacts with this struct */
typedef struct SCS_A_DATA_MATRIX ScsMatrix;
/* stores the necessary private workspace, only the linear system solver
 * interacts with this struct */
typedef struct SCS_LIN_SYS_WORK ScsLinSysWork;

typedef struct SCS_PROBLEM_DATA ScsData;
typedef struct SCS_SETTINGS ScsSettings;
typedef struct SCS_SOL_VARS ScsSolution;
typedef struct SCS_INFO ScsInfo;
typedef struct SCS_SCALING ScsScaling;
typedef struct SCS_WORK ScsWork;
typedef struct SCS_RESIDUALS ScsResiduals;
typedef struct SCS_CONE ScsCone;
typedef struct SCS_ACCEL_WORK ScsAccelWork;
typedef struct SCS_CONE_WORK ScsConeWork;

/* struct containing problem data */
struct SCS_PROBLEM_DATA {
  /* these cannot change for multiple runs for the same call to SCS(init) */
  scs_int m, n; /* A has m rows, n cols */
  ScsMatrix *A; /* A is supplied in data format specified by linsys solver */

  /* these can change for multiple runs for the same call to SCS(init) */
  scs_float *b, *c; /* dense arrays for b (size m), c (size n) */

  ScsSettings *stgs; /* contains solver settings specified by user */
};

/* ScsSettings struct */
struct SCS_SETTINGS {
  /* settings parameters: default suggested input */

  /* these *cannot* change for multiple runs with the same call to SCS(init) */
  scs_int normalize; /* boolean, heuristic data rescaling: 1 */
  scs_float scale;   /* if normalized, rescales by this factor: 5 */
  scs_float rho_x;   /* x equality constraint scaling: 1e-3 */

  /* these can change for multiple runs with the same call to SCS(init) */
  scs_int max_iters;  /* maximum iterations to take: 2500 */
  scs_float eps;      /* convergence tolerance: 1e-3 */
  scs_float alpha;    /* relaxation parameter: 1.8 */
  scs_float cg_rate;  /* for indirect, tolerance goes down like
                         (1/iter)^cg_rate: 2 */
  scs_int verbose;    /* boolean, write out progress: 1 */
  scs_int warm_start; /* boolean, warm start (put initial guess in ScsSolution
                         struct): 0 */
  scs_int acceleration_lookback; /* memory for acceleration */
  const char* write_data_filename; /* string, if set will dump data */
};

/* NB: rows of data matrix A must be specified in this exact order */
struct SCS_CONE {
  scs_int f;     /* number of linear equality constraints */
  scs_int l;     /* length of LP cone */
  scs_int *q;    /* array of second-order cone constraints */
  scs_int qsize; /* length of SOC array */
  scs_int *s;    /* array of SD constraints */
  scs_int ssize; /* length of SD array */
  scs_int ep;    /* number of primal exponential cone triples */
  scs_int ed;    /* number of dual exponential cone triples */
  scs_float *p;  /* array of power cone params, must be \in [-1, 1],
                    negative values are interpreted as specifying the
                    dual cone */
  scs_int psize; /* number of (primal and dual) power cone triples */
};

/* contains primal-dual solution arrays */
struct SCS_SOL_VARS {
  scs_float *x, *y, *s;
};

/* contains terminating information */
struct SCS_INFO {
  scs_int iter;         /* number of iterations taken */
  char status[32];      /* status string, e.g. 'Solved' */
  scs_int status_val;   /* status as scs_int, defined in glbopts.h */
  scs_float pobj;       /* primal objective */
  scs_float dobj;       /* dual objective */
  scs_float res_pri;    /* primal equality residual */
  scs_float res_dual;   /* dual equality residual */
  scs_float res_infeas; /* infeasibility cert residual */
  scs_float res_unbdd;  /* unbounded cert residual */
  scs_float rel_gap;    /* relative duality gap */
  scs_float setup_time; /* time taken for setup phase (milliseconds) */
  scs_float solve_time; /* time taken for solve phase (milliseconds) */
};

/* contains normalization variables */
struct SCS_SCALING {
  scs_float *D, *E; /* for normalization */
  scs_float mean_norm_row_a, mean_norm_col_a;
};

/*
 * main library api's:
 * SCS(init): allocates memory etc (e.g., factorize matrix [I A; A^T -I])
 * SCS(solve): can be called many times with different b,c data per init call
 * SCS(finish): cleans up the memory (one per init call)
 */
ScsWork *SCS(init)(const ScsData *d, const ScsCone *k, ScsInfo *info);
scs_int SCS(solve)(ScsWork *w, const ScsData *d, const ScsCone *k,
                   ScsSolution *sol, ScsInfo *info);
void SCS(finish)(ScsWork *w);
/* scs calls SCS(init), SCS(solve), and SCS(finish) */
scs_int scs(const ScsData *d, const ScsCone *k, ScsSolution *sol,
            ScsInfo *info);

const char *SCS(version)(void);
size_t SCS(sizeof_int)(void);
size_t SCS(sizeof_float)(void);

/* the following structs are not exposed to user */

/* workspace for SCS */
struct SCS_WORK {
  /* x_prev = x from previous iteration */
  scs_float *u, *u_best, *v, *v_best, *u_t, *u_prev, *v_prev;
  scs_float *h, *g, *pr, *dr;
  scs_float g_th, sc_b, sc_c, nm_b, nm_c, best_max_residual;
  scs_float *b, *c;       /* (possibly normalized) b and c vectors */
  scs_int m, n;           /* A has m rows, n cols */
  ScsMatrix *A;           /* (possibly normalized) A matrix */
  ScsLinSysWork *p;       /* struct populated by linear system solver */
  ScsSettings *stgs;      /* contains solver settings specified by user */
  ScsScaling *scal;       /* contains the re-scaling data */
  ScsConeWork *cone_work; /* workspace for the cone projection step */
  AaWork *accel;          /* Struct for acceleration workspace */
};

/* to hold residual information (unnormalized) */
struct SCS_RESIDUALS {
  scs_int last_iter;
  scs_float res_dual;
  scs_float res_pri;
  scs_float res_infeas;
  scs_float res_unbdd;
  scs_float rel_gap;
  scs_float ct_x_by_tau; /* not divided by tau */
  scs_float bt_y_by_tau; /* not divided by tau */
  scs_float tau;
  scs_float kap;
};

#ifdef __cplusplus
}
#endif
#endif
