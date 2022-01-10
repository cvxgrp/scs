#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include <math.h>

void SCS(scale_array)(scs_float *a, const scs_float b, scs_int len);
scs_float SCS(dot)(const scs_float *x, const scs_float *y, scs_int len);
scs_float SCS(norm_sq)(const scs_float *v, scs_int len);
scs_float SCS(norm_2)(const scs_float *v, scs_int len);
scs_float SCS(norm_inf)(const scs_float *a, scs_int l);
void SCS(add_scaled_array)(scs_float *a, const scs_float *b, scs_int n,
                           const scs_float sc);
scs_float SCS(norm_diff)(const scs_float *a, const scs_float *b, scs_int l);
scs_float SCS(norm_inf_diff)(const scs_float *a, const scs_float *b, scs_int l);
scs_float SCS(mean)(const scs_float *x, scs_int l);

#ifdef __cplusplus
}
#endif
#endif
