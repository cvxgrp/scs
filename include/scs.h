#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include <string.h>
#include "cones.h"
#include "linAlg.h"
#include "linSys.h"
#include "util.h"
#include "ctrlc.h"
#include "constants.h"

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
};

/* contains primal-dual solution arrays */
struct SCS_SOL_VARS {
    scs_float *x, *y, *s;
};

/* contains terminating information */
struct SCS_INFO {
    scs_int iter;        /* number of iterations taken */
    char status[32];     /* status string, e.g. 'Solved' */
    scs_int statusVal;   /* status as scs_int, defined in constants.h */
    scs_float pobj;      /* primal objective */
    scs_float dobj;      /* dual objective */
    scs_float resPri;    /* primal equality residual */
    scs_float resDual;   /* dual equality residual */
    scs_float resInfeas; /* infeasibility cert residual */
    scs_float resUnbdd;  /* unbounded cert residual */
    scs_float relGap;    /* relative duality gap */
    scs_float setupTime; /* time taken for setup phase (milliseconds) */
    scs_float solveTime; /* time taken for solve phase (milliseconds) */
};

/* contains normalization variables */
struct SCS_SCALING {
    scs_float *D, *E; /* for normalization */
    scs_float meanNormRowA, meanNormColA;
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
    scs_float *u, *v, *u_t, *u_prev; /* u_prev = u from previous iteration */
    scs_float *h, *g, *pr, *dr;
    scs_float gTh, sc_b, sc_c, nm_b, nm_c;
    scs_float *b, *c;   /* (possibly normalized) b and c vectors */
    scs_int m, n;       /* A has m rows, n cols */
    AMatrix *A;         /* (possibly normalized) A matrix */
    Priv *p;            /* struct populated by linear system solver */
    Settings *stgs;     /* contains solver settings specified by user */
    Scaling *scal;      /* contains the re-scaling data */
    ConeWork *coneWork; /* workspace for the cone projection step */
};

/* to hold residual information (unnormalized) */
struct residuals {
    scs_int lastIter;
    scs_float resDual;
    scs_float resPri;
    scs_float resInfeas;
    scs_float resUnbdd;
    scs_float relGap;
    scs_float cTx_by_tau; /* not divided by tau */
    scs_float bTy_by_tau; /* not divided by tau */
    scs_float tau;
    scs_float kap;
};

#ifdef __cplusplus
}
#endif
#endif
