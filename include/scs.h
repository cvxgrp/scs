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
        AMatrix *A; /* A is supplied in data format specified by linsys solver */

        /* these can change for multiple runs for the same call to scs_init */
        scs_float *b, *c; /* dense arrays for b (size m), c (size n) */

        Settings *stgs; /* contains solver settings specified by user */
    };

    /* Settings struct */
    struct SCS_SETTINGS {
        /* settings parameters: default suggested input */

        /* these *cannot* change for multiple runs with the same call to scs_init */
        scs_int normalize; /* boolean, heuristic data rescaling: 1 */
        scs_float scale; /* if normalized, rescales by this factor: 5 */
        scs_float rho_x; /* x equality constraint scaling: 1e-3 */

        /* these can change for multiple runs with the same call to scs_init */
        scs_int max_iters; /* maximum iterations to take: 2500 */
        scs_float eps; /* convergence tolerance: 1e-3 */
        scs_float alpha; /* relaxation parameter: 1.8 */
        scs_float cg_rate; /* for indirect, tolerance goes down like
                           (1/iter)^cg_rate: 2 */
        scs_int verbose; /* boolean, write out progress: 1 */
        scs_int warm_start; /* boolean, warm start (put initial guess in Sol
                           struct): 0 */
        /* superscs */
        scs_int k0; /* boolean, K0: 1 */
        scs_float c_bl; /* parameter for blind updates: 0.999 */
        scs_int k1; /* boolean, K1: 1 */
        scs_int k2; /* boolean, K2: 1 */
        scs_int nominal; /* boolean, nominal updates: 1 */

        /* line-search */
        scs_int ls; /* max line-search iterations */
        scs_float beta; /* stepsize reduction */
        scs_float sigma; /* line-search parameter */

        /* direction */
        scs_int direction; /* choice of direction: 1 for L-Broyden */
        scs_int tRule; /* rule for selecting relaxation parameter */
        scs_float delta; /* parameter in Broyden realaxation: 0.5 */
        scs_float thetabar; /* modified Broyden's parameter: 1e-1 */
        scs_float alphaC; /* parameter for skipping rule: 1e-2 */
        scs_int memory; /* memory for limited memory QN: 10 */
        scs_int sc_init; /* Boolean, initial scaling for QN: 0 */
    };

    /* contains primal-dual solution arrays */
    struct SCS_SOL_VARS {
        scs_float *x, *y, *s;
    };

    /* contains terminating information */
    struct SCS_INFO {
        scs_int iter; /* number of iterations taken */
        char status[32]; /* status string, e.g. 'Solved' */
        scs_int statusVal; /* status as scs_int, defined in constants.h */
        scs_float pobj; /* primal objective */
        scs_float dobj; /* dual objective */
        scs_float resPri; /* primal equality residual */
        scs_float resDual; /* dual equality residual */
        scs_float resInfeas; /* infeasibility cert residual */
        scs_float resUnbdd; /* unbounded cert residual */
        scs_float relGap; /* relative duality gap */
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


    struct SCS_SU_MEMORY {
        scs_float *S; /**< cached values of s_i (s-memory) */
        scs_float *U; /**< cached values of u_i = (s_i - s_tilde_i)/(s_i'*s_tilde_i) (u-memory)*/
        scs_int mem_current; /**< current memory (before it's full) */
        scs_int mem_idx; /**< head-index of memory */
        scs_int mem; /**< (target) memory */
    };

    /**
     * A finite-memory cache where (Y, S) are cached together with their
     * inner products YS = Y'*S.
     */
    typedef struct SCS_SU_MEMORY SUCache;

    /* workspace for SCS */
    struct SCS_WORK {
        scs_int m; /**< row dimension of A */
        scs_int n; /**< column dimension of A */
        scs_int l; /**< l = m + n + 1 (length of Yk, Sk, etc)*/

        scs_float *u;
        scs_float *v;
        scs_float *u_t;
        scs_float * u_prev;
        scs_float *u_b; /**< u_prev = u from previous iteration */
        scs_float *h;
        scs_float *g;
        scs_float *pr;
        scs_float *dr;
        scs_float gTh;
        scs_float sc_b;
        scs_float sc_c;
        scs_float nm_b;
        scs_float nm_c;
        scs_float *b, *c; /**<  (scpossibly normalized) b and c vectors */
        scs_float *R, *sc_R, *sc_R_prev; /**<  fixed point residuals */
        scs_float *dir, *dut; /**<  variables for direction */
        scs_float *wu, *sc_Rwu; /**< from line search */
        scs_float nrmR_con; /**< \|R\| */
        scs_float *Sk; /**< Sk */
        scs_float *Yk; /**< Yk */

        /**
         *  The (possibly normalized) A matrix 
         */
        AMatrix *A;

        /** 
         * struct populated by linear system solver 
         */
        Priv *p;

        /** 
         * contains solver settings specified by user 
         */
        Settings *stgs;

        /**
         *  contains the re-scaling data 
         */
        Scaling *scal;

        /** 
         * workspace for the cone projection step 
         */
        ConeWork *coneWork;

        /**
         * A cache of Y and S (used, for instance, to compute Broyden-type 
         * or other quasi-Newton directions).
         */
        SUCache *ys_cache;
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
