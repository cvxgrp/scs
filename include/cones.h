#ifndef CONES_H_GUARD
#define CONES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs_blas.h"

/* NB: rows of data matrix A must be specified in this exact order */
struct SCS_CONE {
    scs_int f;     /* number of linear equality constraints */
    scs_int l;     /* length of LP cone */
    scs_int *q;    /* array of second-order cone constraints */
    scs_int qsize; /* length of SOC array */
    scs_int *s;    /* array of SD constraints */
    scs_int ssize; /* length of SD array */
    scs_int ep;    /* number of primal exponential cone triples */
    scs_int ed;    /* number of dual exponential cone triples */
    scs_float *
        p;         /* array of power cone params, must be \in [-1, 1],
                     negative values are interpreted as specifying the dual cone */
    scs_int psize; /* number of (primal and dual) power cone triples */
};

/* private data to help cone projection step */
typedef struct {
#ifdef LAPACK_LIB_FOUND
    /* workspace for eigenvector decompositions: */
    scs_float *Xs, *Z, *e, *work;
    blasint *iwork, lwork, liwork;
#endif
    scs_float totalConeTime;
} ConeWork;

/*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size larger
 * than 1
 * returns length of boundaries array, boundaries malloc-ed here so should be
 * freed
 */
scs_int getConeBoundaries(const Cone *k, scs_int **boundaries);

ConeWork *initCone(const Cone *k);
char *getConeHeader(const Cone *k);
scs_int validateCones(const Data *d, const Cone *k);

/* pass in iter to control how accurate the cone projection
 with iteration, set iter < 0 for exact projection, warm_start contains guess
 of solution, can be SCS_NULL*/
scs_int projDualCone(scs_float *x, const Cone *k, ConeWork *c,
                     const scs_float *warm_start, scs_int iter);
void finishCone(ConeWork *coneWork);
char *getConeSummary(const Info *info, ConeWork *c);

#ifdef __cplusplus
}
#endif
#endif
