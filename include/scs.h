#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#include "glbopts.h"
#include <string.h>
#include "cones.h"
#include "linAlg.h"
#include "linSys.h"
#include "util.h"
#include "ctrlc.h"
#include "constants.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * A cache of \f$(s_i, u_i)\f$ where \f$u_i = \frac{s_i - \tilde{s}_i}{\langle s_i, \tilde{s}_i\rangle}\f$.
     * We do not need to store past values of \f$\tilde{s}_i\f$.
     */
    struct SCS_SU_MEMORY {
        scs_float *S; /**< cached values of \f$s_i\f$ (s-memory) */
        scs_float *U; /**< cached values of \f$u_i = \frac{s_i - \tilde{s}_i}{\langle s_i, \tilde{s}_i\rangle}\f$ (u-memory)*/
        scs_int mem_cursor; /**< current memory cursor [0..mem-1] */
        scs_int mem; /**< (target) memory */
    };    
    
    
    /**
     *  Workspace for SCS 
     */
    struct SCS_WORK {
        scs_int m; /**< row dimension of \f$A\f$ */
        scs_int n; /**< column dimension of \f$A\f$ */
        scs_int l; /**< \f$l = m + n + 1\f$ (length of \f$y_k\f$, \f$s_k\f$, etc)*/
        scs_float *u;
        scs_float *v;
        scs_float *u_t;
        scs_float * u_prev;/**< u_prev = u from previous iteration */
        scs_float *u_b; 
        scs_float *h;
        scs_float *g;
        scs_float *pr;
        scs_float *dr;
        scs_float gTh;
        scs_float sc_b;
        scs_float sc_c;
        scs_float nm_b;
        scs_float nm_c;
        scs_float *b;/**<  (possibly normalized) \c b vector */
        scs_float *c; /**<  (possibly normalized) \c c vector */
        scs_float *R; /**< FPR \f$Ru_k\f$ (fixed point residual) */
        scs_float *R_prev; 
        scs_float *dir; /**< direction \f$d_k\f$ for \f$u\f$ */
        scs_float *dut; /**<  direction \f$\tilde{u} \f$*/
        scs_float *wu; 
        scs_float *wu_t; 
        scs_float *wu_b; 
        scs_float *Rwu; /**< from line search */
        scs_float nrmR_con; /**<  \f$\|Ru_k\|\f$ */
        scs_float *Sk; /**< \f$s_k\f$ */
        scs_float *Yk; /**< \f$y_k\f$ */
        scs_float stepsize; /**< The current stepsize */
        scs_float *s_b; /**< primal residual */
        scs_float kap_b; /**< variable for cerfificates of infeasibility/unboudedness */

        /**
         *  The (possibly normalized) \c A matrix 
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
         * A cache of \c S and \c U (used, to compute Broyden-type directions).
         */
        SUCache *su_cache;
    };

    /**
     *  Structure to hold residual information (unnormalized) 
     */
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
    
    /**
     *  struct containing problem data 
     */
    struct SCS_PROBLEM_DATA {
        /* these cannot change for multiple runs for the same call to scs_init */
        scs_int m, n; /* A has m rows, n cols */
        AMatrix *A; /**< \c A is supplied in data format specified by linsys solver */

        /* these can change for multiple runs for the same call to scs_init */
        scs_float *b, *c; /* dense arrays for b (size m), c (size n) */

        Settings *stgs; /**< contains solver settings specified by user */
    };

    /**
     *  Settings struct 
     */
    struct SCS_SETTINGS {
        /* settings parameters: default suggested input */

        /* these *cannot* change for multiple runs with the same call to scs_init */
        scs_int normalize; /**< boolean, heuristic data rescaling: 1 */
        scs_float scale; /**< if normalized, rescales by this factor: 5 */
        scs_float rho_x; /**< x equality constraint scaling: 1e-3 */

        /* these can change for multiple runs with the same call to scs_init */
        scs_int max_iters; /**< maximum iterations to take: 2500 */
        scs_float eps; /**< convergence tolerance: 1e-3 */
        scs_float alpha; /**< relaxation parameter: 1.8 */
        scs_float cg_rate; /**< for indirect, tolerance goes down like
                           (1/iter)^cg_rate: 2 */
        scs_int verbose; /**< boolean, write out progress: 1 */
        scs_int warm_start; /**< boolean, warm start (put initial guess in Sol
                           struct): 0 */
        /* superscs */
        scs_int do_super_scs; /**< boolean: whether to use superscs or not */
        scs_int k0; /**< boolean, K0: 1 */
        scs_float c_bl; /**< parameter for blind updates: 0.999 */
        scs_int k1; /**< boolean, K1: 1 */
        scs_int k2; /**< boolean, K2: 1 */
        scs_int nominal; /**< boolean, nominal updates: 1 */
        scs_float c1; /**< Parameter to check condition at K1 */
        scs_float sse; /**< Parameter to update r_safe at K1 (denoted as 'q' in the paper) */
        
        /* line-search */
        scs_int ls; /**< max line-search iterations */
        scs_float beta; /**< stepsize reduction */
        scs_float sigma; /**< line-search parameter */

        /* direction */
        direction_type direction; /**< choice of direction: 1 for L-Broyden */
        scs_int tRule; /**< rule for selecting relaxation parameter */
        scs_float delta; /**< parameter in Broyden realaxation: 0.5 */
        scs_float thetabar; /**< modified Broyden's parameter: 1e-1 */
        scs_float alphaC; /**< parameter for skipping rule: 1e-2 */
        scs_int memory; /**< memory for limited memory QN: 10 */
        scs_int sc_init; /**< Boolean, initial scaling for QN: 0 */
    };

    /**
     *  contains primal-dual solution arrays 
     */
    struct SCS_SOL_VARS {
        scs_float *x, *y, *s;
    };

    /**
     *  contains terminating information 
     */
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

    /**
     *  contains normalization variables 
     */
    struct SCS_SCALING {
        scs_float *D, *E; /* for normalization */
        scs_float meanNormRowA, meanNormColA;
    };
    
    /**
     * Creates a new empty solution structure which is to be used 
     * to retrieve the solution \f$(x^\star, y^\star, s^\star)\f$. 
     * 
     * This function does not initialize of allocate memory for 
     * 
     * @return Initialized Sol structure.
     */
    Sol * initSol(void);
    
    /**
     * Creates a new empty #Info structure which is then provided to #scs to get 
     * information about the status of the algorithm (e.g., the duality gap, 
     * the solution status, etc).
     * 
     * @return Initialized #Info structure.
     */
    Info * initInfo(void);
    
    
    /**
     * Creates a new \c #Data structure without allocating memory for \f$A\f$, 
     * \f$b\f$ and \f$c\f$ and its sets all settings to their default values, that 
     * is
     * 
     *  1. alpha = 0.5
     *  2. c0 = 0.999
     *  3. c1 = ...
     * 
     * @return 
     */
    Data * initData(void);
    
    /** 
     * scs calls \c scs_init, \c scs_solve, and \c scs_finish 
     */
    scs_int scs(
            const Data *d, 
            const Cone *k, 
            Sol *sol, 
            Info *info);
    
    /**
     * Returns the version of SCS
     * @return 
     */
    const char *scs_version(void);

#ifdef __cplusplus
}
#endif
#endif
