#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#include "glbopts.h"
#include <string.h>
#include <math.h>
#include "cones.h"
#include "linAlg.h"
#include "linSys.h"
#include "util.h"

/* struct that containing standard problem data */
struct PROBLEM_DATA {
	/* problem dimensions */
	idxint m, n; /* A has m rows, n cols*/

	AMatrix * A; /* A is supplied in data format specified by linsys solver */
	pfloat * b, *c; /* dense arrays for b (size m), c (size n) */

	/* other input parameters: default suggested input */
	idxint MAX_ITERS; /* maximum iterations to take: 2500 */
	pfloat EPS; /* convergence tolerance: 1e-3 */
	pfloat ALPHA; /* relaxation parameter: 1.8 */
	pfloat RHO_X; /* x equality constraint scaling: 1e-3 */
	pfloat SCALE; /* if normalized, rescales by this factor: 5 */
	pfloat CG_RATE; /* for indirect, tolerance goes down like (1/iter)^CG_RATE: 1.5 */
	idxint VERBOSE; /* boolean, write out progress: 1 */
	idxint NORMALIZE; /* boolean, heuristic data rescaling: 1 */
	idxint WARM_START; /* boolean, warm start (put initial guess in Sol struct): 0 */
};

/* contains primal-dual solution arrays */
struct SOL_VARS {
	pfloat * x, *y, *s;
};

/* contains terminating information */
struct INFO {
	idxint iter; /* number of iterations taken */
	char status[32]; /* status string, e.g. 'Solved' */
	idxint statusVal; /* status as idxint, defined below */
	pfloat pobj; /* primal objective */
	pfloat dobj; /* dual objective */
	pfloat resPri; /* primal equality residual */
	pfloat resDual; /* dual equality residual */
	pfloat relGap; /* relative duality gap */
	pfloat setupTime; /* time taken for setup phase */
	pfloat solveTime; /* time taken for solve phase */
};

/* scs returns one of the following integers: (zero should never be returned) */
#define FAILURE -4
#define INDETERMINATE -3
#define INFEASIBLE -2 /* primal infeasible, dual unbounded */
#define UNBOUNDED -1 /* primal unbounded, dual infeasible */
#define SOLVED 1

/* main library api's:
 scs_init: allocates memory (direct version factorizes matrix [I A; A^T -I])
 scs_solve: can be called many times with different b,c data for one init call
 scs_finish: cleans up the memory (one per init call)
 */
Work * scs_init(Data * d, Cone * k, Info * info);
idxint scs_solve(Work * w, Data * d, Cone * k, Sol * sol, Info * info);
void scs_finish(Data * d, Work * w);
/* scs calls scs_init, scs_solve, and scs_finish */
idxint scs(Data * d, Cone * k, Sol * sol, Info * info);

/* the following structs do not need to be exposed */
struct WORK {
	pfloat *u, *v, *u_t, *u_prev;
	pfloat *h, *g, *pr, *dr;
	pfloat gTh, sc_b, sc_c, nm_b, nm_c, meanNormRowA;
	pfloat *D, *E; /* for normalization */
	Priv * p;
};

/* to hold residual information */
struct residuals {
	pfloat resDual;
	pfloat resPri;
	pfloat relGap;
	pfloat cTx;
	pfloat bTy;
	pfloat tau;
	pfloat kap;
};

#endif
