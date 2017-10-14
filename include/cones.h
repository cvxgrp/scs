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
  scs_float *p;  /* array of power cone params, must be \in [-1, 1],
                    negative values are interpreted as specifying the
                    dual cone */
  scs_int psize; /* number of (primal and dual) power cone triples */
};

/* private data to help cone projection step */
typedef struct {
  scs_float total_cone_time;
#ifdef LAPACK_LIB_FOUND
  /* workspace for eigenvector decompositions: */
  scs_float *Xs, *Z, *e, *work;
  blasint *iwork, lwork, liwork;
#endif
} Cone_work;

/*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size larger
 * than 1
 * returns length of boundaries array, boundaries malloc-ed here so should be
 * freed
 */
scs_int get_cone_boundaries(const Cone *k, scs_int **boundaries);

Cone_work *init_cone(const Cone *k);
char *get_cone_header(const Cone *k);
scs_int validate_cones(const Data *d, const Cone *k);

/* pass in iter to control how accurate the cone projection
 with iteration, set iter < 0 for exact projection, warm_start contains guess
 of solution, can be SCS_NULL*/
scs_int proj_dual_cone(scs_float *x, const Cone *k, Cone_work *c,
                       const scs_float *warm_start, scs_int iter);
void finish_cone(Cone_work *c);
char *get_cone_summary(const Info *info, Cone_work *c);

#ifdef __cplusplus
}
#endif
#endif
