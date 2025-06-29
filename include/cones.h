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

#ifdef USE_SPECTRAL_CONES
#include "util_spectral_cones.h" // for newton_stats

// macro for time measurements of SpectralSCS
#ifdef SPECTRAL_TIMING_FLAG
#define SPECTRAL_TIMING(action) action
#else
#define SPECTRAL_TIMING(action)
#endif
#endif

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
  scs_float *Xs, *Z, *e, *work, *rwork;
  scs_complex_float *cXs, *cZ, *cwork;
  blas_int *isuppz, *iwork;
  blas_int lwork, lcwork, lrwork, liwork;
#endif

#ifdef USE_SPECTRAL_CONES
  /* if the projection onto the logarithmic cone should be warmstarted*/
  bool *log_cone_warmstarts;

  /* Needed for ell1 norm cone projection */
  Value_index *work_ell1;
  scs_float *work_ell1_proj;

  // used for timing spectral vector cone and spectral matrix cone projections
  SPECTRAL_TIMING(scs_float tot_time_mat_cone_proj;)
  SPECTRAL_TIMING(scs_float tot_time_vec_cone_proj;)

  /* workspace for singular value decompositions: */
  scs_float *s_nuc, *u_nuc, *vt_nuc, *work_nuc;
  blas_int lwork_nuc;

  /* workspace that is used internally in the logdet projection (for example,
     the gradient and Hessian of the objective function in the projection
     problem are stored using this memory) */
  scs_float *work_logdet;

  /* workspace to store the projection onto the logarithm cone */
  scs_float *saved_log_projs;

  /* Stats for spectral projections, assuming there is only one spectral cone */
  Newton_stats newton_stats;

  /* workspace for projection onto sum-largest-evals cone */
  scs_float *work_sum_of_largest;
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
