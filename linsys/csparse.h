/* Routines modified from CSparse, T. Davis et al */

#ifndef CS_H_GUARD
#define CS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"

ScsMatrix *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                           scs_int triplet);
ScsMatrix *SCS(cs_done)(ScsMatrix *C, void *w, void *x, scs_int ok);
ScsMatrix *SCS(cs_compress)(const ScsMatrix *T, scs_int nz,
                            scs_int *idx_mapping);
ScsMatrix *SCS(cs_spfree)(ScsMatrix *A);
scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);
/* Forms KKT matrix */
ScsMatrix *SCS(form_kkt)(const ScsMatrix *A, const ScsMatrix *P,
                         scs_float *diag_p, const scs_float *diag_r,
                         scs_int *diag_r_idxs, scs_int upper);
#ifdef __cplusplus
}
#endif
#endif
