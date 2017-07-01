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
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * \brief A cache of \f$(s_i, u_i)\f$ where \f$u_i = \frac{s_i - \tilde{s}_i}{\langle s_i, \tilde{s}_i\rangle}\f$.
     * 
     * \note We do not need to store past values of \f$\tilde{s}_i\f$.
     */
    struct SCS_SU_MEMORY {
        scs_float *S; /**< \brief cached values of \f$s_i\f$ (s-memory) */
        scs_float *U; /**< \brief cached values of \f$u_i = \frac{s_i - \tilde{s}_i}{\langle s_i, \tilde{s}_i\rangle}\f$ (u-memory)*/
        scs_int mem_cursor; /**< \brief current memory cursor [0..mem-1] */
        scs_int mem; /**< \brief (target/maximum/allocated) memory */
    };

    /**
     *  \brief Workspace for SCS 
     */
    struct SCS_WORK {
        /**
         *  \brief Row dimension of \f$A\f$. 
         */
        scs_int m;
        /**
         *  \brief Column dimension of \f$A\f$.
         */
        scs_int n;
        /**
         *  \brief Dimension of \f$u_k\f$, that is \f$l = m + n + 1\f$ .
         */
        scs_int l;
        /**
         * \brief Vector \f$u_k=(x_k,y_k,\tau_k)\f$.
         */
        scs_float *u;
        /**
         * \brief Vector \f$v_k = Qu_k\f$ (used only in SCS, not in SuperSCS).
         */
        scs_float *v;
        /**
         * \brief Vector \f$\tilde{u}_k\f$.
         */
        scs_float *u_t;
        /**
         * \brief Vector \f$u_{k-1}\f$ of the previous iteration.
         */
        scs_float * u_prev;
        /**
         * \brief Vector \f$\bar{u}_k\f$.
         */
        scs_float *u_b;
        /**
         */
        scs_float *h;
        /**
         */
        scs_float *g;
        /**
         * \brief Primal residual vector
         * 
         * \f[
         *  \text{pr} = A x + s - b \tau,
         * \f]
         * and in SuperSCS, \f$x\f$ is \f$\bar{x}\f$, \f$s\f$ is \f$\bar{s}\f$
         * and \f$\tau\f$ is \f$\bar{\tau}\f$
         */
        scs_float *pr;
        /**
         * \brief Dual residual vector
         * 
         * \f[
         *  \text{dr} = A'y + c \tau,
         * \f]
         * and in SuperSCS, \f$y\f$ is \f$\bar{y}\f$ and \f$\tau\f$ is \f$\bar{\tau}\f$.
         */
        scs_float *dr;
        /**
         */
        scs_float gTh;
        /**
         * \brief Scaling factor corresponding to \f$b\f$
         */
        scs_float sc_b;
        /**
         * \brief Scaling factor corresponding to \f$c\f$
         */
        scs_float sc_c;
        /**
         * \brief Norm of \f$b\f$
         */
        scs_float nm_b;
        /**
         * \brief Norm of \f$c\f$
         */
        scs_float nm_c;
        /**  
         * \brief The (possibly normalized) vector \f$b\f$.
         */
        scs_float *b;
        /**  
         * \brief The (possibly normalized) vector \f$c\f$.
         */
        scs_float *c;
        /**  
         * \brief Fixed-point residual \f$R_k\f$.
         */
        scs_float *R;
        /**  
         * \brief Fixed-point residual (FPR) of the previous iteration \f$R_{k-1}\f$.
         */
        scs_float *R_prev;
        /**
         * \brief Direction \f$d_k\f$
         */
        scs_float *dir;
        /**
         * \brief Hessian approximation used by the full Broyden method.
         * 
         * @see ::full_broyden
         */
        scs_float * H;
        /** 
         * \brief Direction corresponding to \f$\tilde{u}\f$.
         */
        scs_float *dut;
        /**
         */
        scs_float *wu;
        /**
         */
        scs_float *wu_t;
        /**
         */
        scs_float *wu_b;
        /**
         * \brief Vector \c Rwu from line search.
         */
        scs_float *Rwu;
        /** 
         * \brief \f$\|Ru_k\|\f$. 
         */
        scs_float nrmR_con;
        /**
         *  \brief \f$s_k = u_k - u_{k-1}\f$ 
         */
        scs_float *Sk;
        /** 
         * \brief \f$y_k = R_k - R_{k-1}\f$ 
         */
        scs_float *Yk;
        /** 
         * \brief The current stepsize \f$t_k\f$
         */
        scs_float stepsize;
        /** 
         * \brief Variable that corresponds to the primal slack for the 2nd step of DRS 
         */
        scs_float *s_b;
        /** 
         * \brief Variable for certificates of infeasibility/unboudedness 
         */
        scs_float kap_b;
        /**
         *  \brief The (possibly normalized) \c A matrix 
         */
        AMatrix *A;
        /** 
         * \brief struct populated by linear system solver 
         */
        Priv *p;
        /** 
         * \brief contains solver settings specified by user 
         */
        Settings *stgs;
        /**
         * \brief contains the re-scaling data 
         */
        Scaling *scal;
        /** 
         * \brief workspace for the cone projection step 
         */
        ConeWork *coneWork;
        /**
         * \brief A cache of \c S and \c U (used, to compute Broyden-type directions).
         */
        SUCache *su_cache;
    };

    /**
     *  \brief Structure to hold residual information (unnormalized) 
     */
    struct residuals {
        /**
         * \brief The last iteration when the residuals were updated.
         */
        scs_int lastIter;
        /**
         * \brief Dual residual
         * 
         * \f[
         * \text{resdual} = \frac{E(A'y + \tau c)}{\tau(1+\|c\|)\cdot \text{scale}_c\cdot \text{scale}}
         * \f]
         */
        scs_float resDual;
        /**
         * \brief Primal residual
         * 
         * \f[
         *  \text{respri} = \frac{\|D(Ax+s-\tau b)\|}{\tau(1+\|b\|)(\text{scale}_b\cdot \text{scale})}
         * \f]
         */
        scs_float resPri;
        /**
         * \brief Infeasibility residual
         * 
         * \f[
         *  \text{infres} = -\frac{\|Db\| \|EA'y\|}{b'y \cdot \text{scale}}
         * \f]
         */
        scs_float resInfeas;
        /**
         * \brief Unboundedness
         * 
         * \f[
         * \text{unbdd} = -\frac{\|Ec\| \|D(Ax+s)}{c'x\cdot \text{scale}}
         * \f]
         */
        scs_float resUnbdd;
        /**
         * \brief Relative duality gap defined as 
         * 
         * \f[
         *  \text{relgap} = \frac{c'x + b'y}{1+|c'x|+|b'y|}
         * \f]
         */
        scs_float relGap;
        scs_float cTx_by_tau; /* not divided by tau */
        scs_float bTy_by_tau; /* not divided by tau */
        /**
         * Variable \f$\tau\f$ (\f$\bar{\tau}\f$ in SuperSCS)
         */
        scs_float tau; /* for superSCS it's tau_b */
        /**
         * Variable \f$\kappa\f$ (\f$\bar{\kappa}\f$ in SuperSCS)
         */
        scs_float kap; /* for superSCS it's kap_b */
    };

    /**
     *  \brief struct containing problem data 
     */
    struct SCS_PROBLEM_DATA {
        /* these cannot change for multiple runs for the same call to scs_init */
        /**
         *  row dimension of \f$A\f$ 
         */
        scs_int m;
        /**
         *  column dimension of \f$A\f$
         */
        scs_int n;
        AMatrix *A; /**< \c A is supplied in data format specified by linsys solver */

        /* these can change for multiple runs for the same call to scs_init */
        scs_float *b, *c; /* dense arrays for b (size m), c (size n) */

        Settings *stgs; /**< contains solver settings specified by user */
    };

    /**
     * \brief Settings structure
     */
    struct SCS_SETTINGS {
        /* settings parameters: default suggested input */


        /* -------------------------------------
         * General Settings
         * 
         * these *cannot* change for multiple runs 
         * with the same call to scs_init
         * ------------------------------------- */

        /** 
         * Boolean, heuristic data rescaling
         * 
         * Default : 1
         */
        scs_int normalize;

        scs_float scale; /**< if normalized, rescales by this factor: 5 */
        scs_float rho_x; /**< x equality constraint scaling: 1e-3 */


        /* -------------------------------------
         * General Settings
         * 
         * these can change for multiple runs with 
         * the same call to scs_init
         * ------------------------------------- */

        /**
         * Maximum iterations to take: 2500 
         */
        scs_int max_iters;
        /**
         * Maximum iterations of the previous invocation to scs.
         * 
         * Used to avoid memory leaks when recording the progress of the algorithm.
         */
        scs_int previous_max_iters;
        /** 
         * Convergence tolerance.
         * 
         * Default: 1e-3 
         */
        scs_float eps;
        /** 
         * Relaxation parameter.
         * 
         * Default : 1.8 
         */
        scs_float alpha;
        scs_float cg_rate; /**< for indirect, tolerance goes down like
                           (1/iter)^cg_rate: 2 */

        /** 
         * Level of verbosity.
         * 
         * Three levels are supported: 0, 1 and 2.
         * 
         * Default : 1
         * 
         */
        scs_int verbose;
        /** 
         * Boolean, warm start (put initial guess in Sol struct): 0 
         */
        scs_int warm_start;

        /* -------------------------------------
         * Settings associated with SuperSCS
         * ------------------------------------- */

        scs_int do_super_scs; /**< boolean: whether to use superscs or not */
        scs_int k0; /**< boolean, K0: 1 */
        scs_float c_bl; /**< parameter for blind updates: 0.999 */
        scs_int k1; /**< boolean, K1: 1 */
        scs_int k2; /**< boolean, K2: 1 */
        scs_float c1; /**< Parameter to check condition at K1 */
        scs_float sse; /**< Parameter to update r_safe at K1 (denoted as \f$q\f$ in the paper) */

        /* -------------------------------------
         * Settings associated with the line 
         * search
         * ------------------------------------- */
        /** 
         * max line-search iterations 
         */
        scs_int ls;
        /**
         * Step size reduction coefficient. 
         * 
         * In every line search iteration, the step size is reduced as 
         * \f$t \leftarrow \beta t\f$.
         */
        scs_float beta;
        /** 
         * Line-search parameter 
         */
        scs_float sigma;

        /* -------------------------------------
         * Settings associated with the direction
         * ------------------------------------- */
        /** 
         * Choice of direction
         * 
         * Default : ::restarted_broyden 
         */
        direction_type direction;
        /** 
         * Modified Broyden parameter.
         * 
         * Default : 1e-1 
         */
        scs_float thetabar;
        /** 
         * Memory length for limited-memory Broyden method
         * 
         * Default: 10 
         */
        scs_int memory;
        /**
         * Option for the Broyden direction.
         * 
         * Default: 3
         * 
         * @see ::computeDirection
         */
        scs_int tRule;
        /**
         * Boolean; whether an initial scaling is desired 
         * in the full Broyden method
         * 
         * Default : 0
         */
        scs_int broyden_init_scaling;
        /**
         * Whether to record progress data when running SuperSCS.
         */
        scs_int do_record_progress;
    };

    /**
     *  \brief Primal-dual solution arrays 
     */
    struct SCS_SOL_VARS {
        scs_float *x;
        scs_float *y;
        scs_float *s;
    };

    /**
     *  \brief Terminating information 
     * 
     * \see ::freeInfo
     */
    struct SCS_INFO {
        scs_int iter; /**< \brief number of iterations taken */
        char status[32]; /**< \brief status string, e.g. 'Solved' */
        scs_int statusVal; /**< \brief status as scs_int, defined in constants.h */
        scs_float pobj; /**< \brief primal objective */
        scs_float dobj; /**< \brief dual objective */
        scs_float resPri; /**< \brief primal equality residual */
        scs_float resDual; /**< \brief dual equality residual */
        scs_float resInfeas; /**< \brief infeasibility cert residual */
        scs_float resUnbdd; /**< \brief unbounded cert residual */
        scs_float relGap; /**< \brief relative duality gap */
        scs_float setupTime; /**< \brief time taken for setup phase (milliseconds) */
        scs_float solveTime; /**< \brief time taken for solve phase (milliseconds) */

        scs_int history_length; /**< \brief how many history entries */
        scs_int * progress_iter; /**< \brief iterations when residulas are recorded */
        scs_float * progress_relgap; /**< \brief relative gap history */
        scs_float * progress_respri; /**< \brief primal residual history */
        scs_float * progress_resdual; /**< \brief dual residual history */
        scs_float * progress_pcost; /**< \brief scaled primal cost history */
        scs_float * progress_dcost; /**< \brief sclaed dual cost history */
        scs_float * progress_norm_fpr; /**< \brief FPR history */
    };

    /**
     *  \brief Normalization variables 
     */
    struct SCS_SCALING {
        scs_float *D, *E; /* for normalization */
        scs_float meanNormRowA, meanNormColA;
    };

    /**
     * Creates a new empty solution structure which is to be used 
     * to retrieve the solution \f$(x^\star, y^\star, s^\star)\f$. 
     * 
     * This function does not initialize of allocate memory for \c x, \c s
     * or \c y (but it sets the respective pointers to ::SCS_NULL).
     * 
     * @return Initialized ::Sol structure.
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
     * @return ::Data object
     * 
     * @see ::setDefaultSettings
     */
    Data * initData(void);

    /** 
     * scs calls \c scs_init, \c scs_solve, and \c scs_finish 
     * 
     * @param d
     * @param k
     * @param sol
     * @param info
     * 
     * @return status code
     * 
     * \remark It is very important that <code>info</code> is created using 
     * ::initInfo.
     */
    scs_int scs(
            const Data *d,
            const Cone *k,
            Sol *sol,
            Info *info);

    /**
     * Returns the version of SCS
     * 
     * @return 
     */
    const char *scs_version(void);

#ifdef __cplusplus
}
#endif
#endif
