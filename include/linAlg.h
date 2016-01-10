#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include <math.h>

void setAsScaledArray(scs_float *x, const scs_float *a, const scs_float b,
                      scs_int len);
void scaleArray(scs_float *a, const scs_float b, scs_int len);
scs_float innerProd(const scs_float *x, const scs_float *y, scs_int len);
scs_float calcNormSq(const scs_float *v, scs_int len);
scs_float calcNorm(const scs_float *v, scs_int len);
scs_float calcNormInf(const scs_float *a, scs_int l);
void addScaledArray(scs_float *a, const scs_float *b, scs_int n,
                    const scs_float sc);
scs_float calcNormDiff(const scs_float *a, const scs_float *b, scs_int l);
scs_float calcNormInfDiff(const scs_float *a, const scs_float *b, scs_int l);

#ifdef __cplusplus
}
#endif
#endif
