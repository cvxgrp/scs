/*
 * Sparse matrix utilities adapted from CSparse (T. Davis).
 * Provides allocation, triplet-to-CSC compression, and KKT matrix formation.
 */

#ifndef CS_H_GUARD
#define CS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

ScsMatrix *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                           scs_int triplet);
ScsMatrix *SCS(cs_done)(ScsMatrix *C, void *w, void *x, scs_int ok);
ScsMatrix *SCS(cs_compress)(const ScsMatrix *T, scs_int nz,
                            scs_int *idx_mapping);
ScsMatrix *SCS(cs_spfree)(ScsMatrix *A);
scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);
/* Forms KKT matrix */
/* (research) SOC block-metric registry: when set, supporting direct
 * solvers form dense upper-triangular SOC blocks in the KKT (2,2) region
 * instead of pure diagonal, and refresh their values from `vals` on every
 * diag_r update. vals is upper-triangle column-major per block, length
 * sum_b q_b(q_b+1)/2, holding the final KKT values (i.e. -W entries).
 * Not thread-safe; single-solver research use only. */
typedef struct {
  scs_int n;
  const scs_int *starts; /* row offsets within the y block */
  const scs_int *sizes;
  const scs_float *vals;
  scs_int refine; /* iterative-refinement request (set by the core when a
                   * boost is active or residuals near tolerance; sticky
                   * so the effective solve operator changes only once) */
} ScsSocMetric;
void SCS(set_soc_metric)(const ScsSocMetric *sm);
const ScsSocMetric *SCS(get_soc_metric)(void);

ScsMatrix *SCS(form_kkt)(const ScsMatrix *A, const ScsMatrix *P,
                         scs_float *diag_p, const scs_float *diag_r,
                         scs_int *diag_r_idxs, scs_int *soc_idxs,
                         scs_int upper);
#ifdef __cplusplus
}
#endif
#endif
