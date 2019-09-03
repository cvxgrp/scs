#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "scs.h"
#include "util.h"

/* this struct defines the data matrix A */
struct SCS_A_DATA_MATRIX {
  /* A is supplied in column compressed format */
  scs_float *x; /* A values, size: NNZ A */
  scs_int *i;   /* A row index, size: NNZ A */
  scs_int *p;   /* A column pointer, size: n+1 */
  scs_int m, n; /* m rows, n cols */
};

void SCS(_accum_by_atrans)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                           const scs_float *x, scs_float *y);
void SCS(_accum_by_a)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                      const scs_float *x, scs_float *y);
void SCS(_normalize_a)(ScsMatrix *A, const ScsSettings *stgs, const ScsCone *k,
                       ScsScaling *scal);
void SCS(_un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                          const ScsScaling *scal);
scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);

#ifdef __cplusplus
}
#endif
#endif
