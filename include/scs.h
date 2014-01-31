#ifndef SCS_H_GUARD 
#define SCS_H_GUARD

#include <string.h>    
#include <sys/time.h>
#include <math.h>
#include "glbopts.h"
#include "cones.h"
#include "linAlg.h"
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
  idxint VERBOSE, NORMALIZE;  /* boolean */
};

/* contains primal-dual solution vectors */
struct SOL_VARS {
  pfloat * x, * y, *s;
};
    
/* contains terminating information */
struct INFO {
	idxint iter;
	char status[16];
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
  pfloat *h, *g;  
  pfloat gTh, sc_b, sc_c, scale;
  pfloat nm_b, nm_c, nm_Q;
  pfloat *D, *E;
  idxint l;
  char * method;
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

/* these are actually library "api"'s */
idxint scs(Data * d, Cone * k, Sol * sol, Info * info);
void freeData(Data *d, Cone *k);
void freeSol(Sol *sol);
void printSol(Data * d, Sol * sol, Info * info);

#endif
