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
 * solvers carry each block's -W = -rI + 2r e0 e0' - 2r ww' as the scalar
 * diagonal plus two bordered columns appended to the KKT (sparse, exact,
 * quasi-definite; see form_kkt), refreshing border values from `vals` on
 * every diag_r update. vals layout per block, in emission order:
 * [d1 = -1/(2r), w_0 .. w_{q-1}, d2 = +1/(2r)], length sum_b (q_b + 2).
 * Not thread-safe; single-solver research use only. */
typedef struct {
  scs_int n;
  const scs_int *starts; /* row offsets within the y block */
  const scs_int *sizes;
  const scs_float *vals;
} ScsSocMetric;
void SCS(set_soc_metric)(const ScsSocMetric *sm);
const ScsSocMetric *SCS(get_soc_metric)(void);

/* borders: 0 = suppress border columns (large SOC blocks fall back to the
 * scalar diagonal; used to compute the AMD ordering on a pattern that
 * includes the small dense blocks but not the borders), 1 = full hybrid */
ScsMatrix *SCS(form_kkt)(const ScsMatrix *A, const ScsMatrix *P,
                         scs_float *diag_p, const scs_float *diag_r,
                         scs_int *diag_r_idxs, scs_int *soc_idxs,
                         scs_int borders, scs_int upper);
#ifdef __cplusplus
}
#endif
#endif
