#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#include "glbopts.h"
#include <string.h>
#include <math.h>
#include "cones.h"
#include "linAlg.h"
#include "linSys.h"
#include "util.h"


/* SCS VERSION NUMBER ----------------------------------------------    */
#define SCS_VERSION ("1.0.6")

/* SCS returns one of the following integers: (zero never returned)     */
#define FAILURE         (-4)
#define INDETERMINATE   (-3)
#define INFEASIBLE      (-2)    /* primal infeasible, dual unbounded    */
#define UNBOUNDED       (-1)    /* primal unbounded, dual infeasible    */
#define SOLVED          (1)

/* DEFAULT SOLVER PARAMETERS AND SETTINGS --------------------------    */
#define MAX_ITERS       (2500)
#define EPS             (1E-3)
#define ALPHA           (1.5)
#define RHO_X           (1E-3)
#define SCALE           (5.0)
#define CG_RATE         (2.0)
#define VERBOSE         (1)
#define NORMALIZE       (1)
#define WARM_START      (0)

/* struct that containing standard problem data */
struct PROBLEM_DATA {
	/* problem dimensions */
	scs_int m, n; /* A has m rows, n cols*/

	AMatrix * A; /* A is supplied in data format specified by linsys solver */
	scs_float * b, *c; /* dense arrays for b (size m), c (size n) */

	/* other input parameters: default suggested input */
	scs_int max_iters; /* maximum iterations to take: 2500 */
	scs_float eps; /* convergence tolerance: 1e-3 */
	scs_float alpha; /* relaxation parameter: 1.8 */
	scs_float rho_x; /* x equality constraint scaling: 1e-3 */
	scs_float scale; /* if normalized, rescales by this factor: 5 */
	scs_float cg_rate; /* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */
	scs_int verbose; /* boolean, write out progress: 1 */
	scs_int normalize; /* boolean, heuristic data rescaling: 1 */
	scs_int warm_start; /* boolean, warm start (put initial guess in Sol struct): 0 */
};

/* contains primal-dual solution arrays */
struct SOL_VARS {
	scs_float * x, *y, *s;
};

/* contains terminating information */
struct INFO {
	scs_int iter; /* number of iterations taken */
	char status[32]; /* status string, e.g. 'Solved' */
	scs_int statusVal; /* status as scs_int, defined below */
	scs_float pobj; /* primal objective */
	scs_float dobj; /* dual objective */
	scs_float resPri; /* primal equality residual */
	scs_float resDual; /* dual equality residual */
	scs_float relGap; /* relative duality gap */
	scs_float setupTime; /* time taken for setup phase */
	scs_float solveTime; /* time taken for solve phase */
};

/* main library api's:
 scs_init: allocates memory (direct version factorizes matrix [I A; A^T -I])
 scs_solve: can be called many times with different b,c data for one init call
 scs_finish: cleans up the memory (one per init call)
 */
Work * scs_init(Data * d, Cone * k, Info * info);
scs_int scs_solve(Work * w, Data * d, Cone * k, Sol * sol, Info * info);
void scs_finish(Data * d, Work * w);
/* scs calls scs_init, scs_solve, and scs_finish */
scs_int scs(Data * d, Cone * k, Sol * sol, Info * info);

/* the following structs do not need to be exposed */
struct WORK {
	scs_float *u, *v, *u_t, *u_prev; /* u_prev = u from previous iteration */
	scs_float *h, *g, *pr, *dr;
	scs_float gTh, sc_b, sc_c, nm_b, nm_c, meanNormRowA, meanNormColA;
	scs_float *D, *E; /* for normalization */
	Priv * p;
};

/* to hold residual information */
struct residuals {
	scs_float resDual;
	scs_float resPri;
	scs_float relGap;
	scs_float cTx;
	scs_float bTy;
	scs_float tau;
	scs_float kap;
};

#endif
