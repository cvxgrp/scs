#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

/* YOUR LINEAR SYSTEM SOLVER MUST IMPLEMENT THESE METHODS AND PRIVATE_DATA STRUCT */

/* private data structs (that you define) containing any necessary data to solve linear system, etc. */
/* this defines the matrix A, only the linear system solver interacts with this struct */
typedef struct A_DATA_MATRIX AMatrix;
/* stores the necessary private workspace, only the linear system solver interacts with this struct */
typedef struct PRIVATE_DATA Priv;

/* initialize Priv structure and perform any necessary preprocessing */
Priv * initPriv(Data * d);
/* solves [d->RHO_X * I  A' ; A  -I] x = b for x, stores result in b, s contains warm-start, iter is current scs iteration count */
scs_int solveLinSys(Data * d, Priv * p, scs_float * b, const scs_float * s, scs_int iter);
/* frees Priv structure and allocated memory in Priv */
void freePriv(Priv * p);

/* forms y += A'*x */
void accumByAtrans(Data * d, Priv * p, const scs_float *x, scs_float *y);
/* forms y += A*x */
void accumByA(Data * d, Priv * p, const scs_float *x, scs_float *y);

/* returns negative num if input data is invalid */
scs_int validateLinSys(Data *d);

/* returns string describing method, can return null, if not null free will be called on output */
char * getLinSysMethod(Data * d, Priv * p);
/* returns string containing summary information about linear system solves, can return null, if not null free will be called on output */
char * getLinSysSummary(Priv * p, Info * info);

/* Normalization routines, used if d->NORMALIZE is true */
/* normalizes A matrix, sets w->E and w->D diagonal scaling matrices, Anew = d->SCALE * (D^-1)*A*(E^-1) (different to paper which is D*A*E)
 * D and E must be all positive entries, D must satisfy cone boundaries
 * must set (w->meanNormRowA = mean of norms of rows of normalized A) THEN scale resulting A by d->SCALE */
void normalizeA(Data * d, Work * w, Cone * k);
/* unnormalizes A matrix, unnormalizes by w->D and w->E and d->SCALE */
void unNormalizeA(Data *d, Work * w);

#endif
