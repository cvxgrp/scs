#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "amatrix.h"
#include "linalg.h"
#include "linsys.h"
#include "scs.h"
#include "util.h"

void SCS(_accum_by_atrans)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                           const scs_float *x, scs_float *y);
void SCS(_accum_by_a)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                      const scs_float *x, scs_float *y);
void SCS(_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                       const ScsCone *k, ScsScaling *scal);
void SCS(_un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                          const ScsScaling *scal);
scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);

#ifdef __cplusplus
}
#endif

#endif
