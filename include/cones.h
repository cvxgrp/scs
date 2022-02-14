#ifndef CONES_H_GUARD
#define CONES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"
#include "scs_blas.h"
#include "scs_work.h"
#include <string.h>

/* private data to help cone projection step */
struct SCS_CONE_WORK {
  /*
   * cone_boundaries will contain array of indices of rows of A corresponding to
   * cone boundaries, boundaries[0] is starting index for cones of size larger
   * than 1
   */
  ScsCone *k; /* original cone information */
  scs_int *cone_boundaries;
  scs_int cone_boundaries_len;
  scs_int scaled_cones; /* boolean, whether the cones have been scaled */
  scs_float *s;         /* used for Moreau decomposition in projection */
  scs_int m;            /* total length of cone */
  /* box cone quantities */
  scs_float box_t_warm_start;
#ifdef USE_LAPACK
  /* workspace for eigenvector decompositions: */
  scs_float *Xs, *Z, *e, *work;
  blas_int lwork;
#endif
};

void SCS(free_cone)(ScsCone *k);
void SCS(deep_copy_cone)(ScsCone *dest, const ScsCone *src);
ScsConeWork *SCS(init_cone)(ScsCone *k, scs_int m);
char *SCS(get_cone_header)(const ScsCone *k);
scs_int SCS(validate_cones)(const ScsData *d, const ScsCone *k);
scs_int SCS(proj_dual_cone)(scs_float *x, ScsConeWork *c, ScsScaling *scal,
                            scs_float *r_y);
void SCS(finish_cone)(ScsConeWork *c);
void SCS(set_r_y)(const ScsConeWork *c, scs_float scale, scs_float *r_y);
void SCS(enforce_cone_boundaries)(const ScsConeWork *c, scs_float *vec,
                                  scs_float (*f)(const scs_float *, scs_int));

#ifdef __cplusplus
}
#endif
#endif
