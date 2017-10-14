#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include <string.h>
#include "cones.h"
#include "lin_alg.h"
#include "lin_sys.h"
#include "util.h"
#include "ctrlc.h"
#include "constants.h"
#include "accel.h"

/* struct containing problem data */
struct SCS_PROBLEM_DATA {
    /* these cannot change for multiple runs for the same call to scs_init */
    scs_int m, n; /* A has m rows, n cols */
    AMatrix *A;   /* A is supplied in data format specified by linsys solver */

    /* these can change for multiple runs for the same call to scs_init */
    scs_float *b, *c; /* dense arrays for b (size m), c (size n) */

    Settings *stgs; /* contains solver settings specified by user */
};

/* Settings struct */
struct SCS_SETTINGS {
    /* settings parameters: default suggested input */

    /* these *cannot* change for multiple runs with the same call to scs_init */
    scs_int normalize; /* boolean, heuristic data rescaling: 1 */
    scs_float scale;   /* if normalized, rescales by this factor: 5 */
    scs_float rho_x;   /* x equality constraint scaling: 1e-3 */

    /* these can change for multiple runs with the same call to scs_init */
    scs_int max_iters;  /* maximum iterations to take: 2500 */
    scs_float eps;      /* convergence tolerance: 1e-3 */
    scs_float alpha;    /* relaxation parameter: 1.8 */
    scs_float cg_rate;  /* for indirect, tolerance goes down like
                           (1/iter)^cg_rate: 2 */
    scs_int verbose;    /* boolean, write out progress: 1 */
    scs_int warm_start; /* boolean, warm start (put initial guess in Sol
                           struct): 0 */
    scs_int acceleration_lookback;
};

/* contains primal-dual solution arrays */
struct SCS_SOL_VARS {
    scs_float *x, *y, *s;
};

/* contains terminating information */
struct SCS_INFO {
    scs_int iter;        /* number of iterations taken */
    char status[32];     /* status string, e.g. 'Solved' */
    scs_int status_val;   /* status as scs_int, defined in constants.h */
    scs_float pobj;      /* primal objective */
    scs_float dobj;      /* dual objective */
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
 * scs_init: allocates memory etc (direct version factorizes matrix [I A; A^T
 * -I])
 * scs_solve: can be called many times with different b,c data for one init call
 * scs_finish: cleans up the memory (one per init call)
 */
Work *scs_init(const Data *d, const Cone *k, Info *info);
scs_int scs_solve(Work *w, const Data *d, const Cone *k, Sol *sol, Info *info);
void scs_finish(Work *w);
/* scs calls scs_init, scs_solve, and scs_finish */
scs_int scs(const Data *d, const Cone *k, Sol *sol, Info *info);
const char *scs_version(void);

/* the following structs are not exposed to user */

/* workspace for SCS */
struct SCS_WORK {
    scs_float *u, *v, *u_t, *u_prev, *v_prev; /* u_prev = u from previous iteration */
    scs_float *h, *g, *pr, *dr;
    scs_float g_th, sc_b, sc_c, nm_b, nm_c;
    scs_float *b, *c;   /* (possibly normalized) b and c vectors */
    scs_int m, n;       /* A has m rows, n cols */
    AMatrix *A;         /* (possibly normalized) A matrix */
    Priv *p;            /* struct populated by linear system solver */
    Accel *accel;       /* Struct for acceleration workspace */
    Settings *stgs;     /* contains solver settings specified by user */
    Scaling *scal;      /* contains the re-scaling data */
    Cone_work *cone_work; /* workspace for the cone projection step */
};

/* to hold residual information (unnormalized) */
struct residuals {
    scs_int last_iter;
    scs_float res_dual;
    scs_float res_pri;
    scs_float res_infeas;
    scs_float res_unbdd;
    scs_float rel_gap;
    scs_float c_tx_by_tau; /* not divided by tau */
    scs_float b_ty_by_tau; /* not divided by tau */
    scs_float tau;
    scs_float kap;
};

#ifdef __cplusplus
}
#endif
#endif
