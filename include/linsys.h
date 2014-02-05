#ifndef LINSYS_H_GUARD 
#define LINSYS_H_GUARD

/* YOUR LINEAR SYSTEM SOLVER MUST IMPLEMENT THESE METHODS AND STRUCTS */

/* private data struct (that you define) containing any necessary data to solve linear system, etc. */
typedef struct PRIVATE_DATA Priv;

/* initialize Priv strucutre and perform any necessary preprocessing */
Priv * initPriv(Data * d);
/* solves [I A';A -I] x = b for x, stores result in b, s contains warm-start, iter is current scs iteration count */
void solveLinSys(Data * d, Priv * p, pfloat * b, const pfloat * s, idxint iter);
/* frees Priv structure and allocated memory in Priv */
void freePriv(Priv * p);

/* forms y += A'*x */
void accumByAtrans(Data * d, Priv * p, const pfloat *x, pfloat *y);
/* forms y += A*x */
void accumByA(Data * d, Priv * p, const pfloat *x, pfloat *y);

/* returns string describing method, can return null, if not null free will be called on output */
char * getLinSysMethod();
/* returns string containing summary information about linear system solves, can return null, if not null free will be called on output */
char * getLinSysSummary(Info * info);

#endif
