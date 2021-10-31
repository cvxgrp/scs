#ifndef CONES_H_GUARD
#define CONES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"
#include "scs_blas.h"

/* private data to help cone projection step */
struct SCS_CONE_WORK {
  /*
   * cone_boundaries will contain array of indices of rows of A corresponding to
   * cone boundaries, boundaries[0] is starting index for cones of size larger
   * than 1
   */
  scs_int *cone_boundaries;
  scs_int cone_boundaries_len;
  scs_int scaled_cones;
  scs_float *s;     /* used for Moreau decomposition in projection */
  scs_int cone_len; /* total length of cone (= m) */
  /* box cone quantities */
  scs_float *bl, *bu, box_t_warm_start;
#ifdef USE_LAPACK
  /* workspace for eigenvector decompositions: */
  scs_float *Xs, *Z, *e, *work;
  blas_int lwork;
#endif
};

ScsConeWork *SCS(init_cone)(const ScsCone *k, scs_int cone_len);
char *SCS(get_cone_header)(const ScsCone *k);
scs_int SCS(validate_cones)(const ScsData *d, const ScsCone *k);
scs_int SCS(proj_dual_cone)(scs_float *x, const ScsCone *k, ScsConeWork *c,
                            ScsScaling *scal, scs_float *rho_y_vec);
void SCS(finish_cone)(ScsConeWork *c);
void SCS(set_rho_y_vec)(const ScsCone *k, const ScsConeWork *c, scs_float scale,
                        scs_float *rho_y_vec);
void SCS(enforce_cone_boundaries)(const ScsCone *k, const ScsConeWork *c,
                                  scs_float *vec);

#ifdef __cplusplus
}
#endif
#endif
