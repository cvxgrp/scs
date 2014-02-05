#ifndef SCS_H_GUARD 
#define SCS_H_GUARD

#include <string.h>    
#include <sys/time.h>
#include <math.h>
#include "glbopts.h"
#include "cones.h"
#include "linAlg.h"
#include "linsys.h"
#include "util.h"

/* struct that containing standard problem data */
struct PROBLEM_DATA {
  idxint n, m; /* problem dimensions */
  /* problem data, A, b, c: */
  pfloat * Ax;
  idxint * Ai, * Ap;
  pfloat * b, * c;
  idxint MAX_ITERS;
  pfloat EPS, ALPHA, UNDET_TOL, RHO_X;
  idxint VERBOSE, NORMALIZE, WARM_START;  /* boolean */
};

/* contains primal-dual solution vectors */
struct SOL_VARS {
  pfloat * x, * y, *s;
};
    
/* contains terminating information */
struct INFO {
	idxint iter;
	char status[32];
	idxint statusVal; /* status as idxint */
    pfloat pobj;
	pfloat dobj;
	pfloat resPri;
	pfloat resDual;
	pfloat relGap;
	pfloat time;
};

/* the following structs do not need to be exposed */
struct WORK {
  pfloat *u, *v, *u_t, *u_prev;
  pfloat *h, *g, *pr, *dr; 
  pfloat gTh, sc_b, sc_c, scale, nm_b, nm_c, meanNormRowA;
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

/* main library api's:
scs_init: allocates memory (direct version factorizes matrix [I A; A^T -I])
scs_solve: can be called many times with different b,c data for one init call
scs_finish: cleans up the memory (one per init call)
*/
Work * scs_init(Data * d, Cone * k);
idxint scs_solve(Work * w, Data * d, Cone * k, Sol * sol, Info * info);
void scs_finish(Data * d, Work * w);
/* scs calls scs_init, scs_solve, and scs_finish */
idxint scs(Data * d, Cone * k, Sol * sol, Info * info);

#endif
