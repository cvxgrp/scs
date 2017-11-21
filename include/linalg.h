#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include <math.h>

void set_as_scaled_array(scs_float *x, const scs_float *a, const scs_float b,
                         scs_int len);
void scale_array(scs_float *a, const scs_float b, scs_int len);
scs_float inner_prod(const scs_float *x, const scs_float *y, scs_int len);
scs_float calc_norm_sq(const scs_float *v, scs_int len);
scs_float calc_norm(const scs_float *v, scs_int len);
scs_float calc_norm_inf(const scs_float *a, scs_int l);
void add_scaled_array(scs_float *a, const scs_float *b, scs_int n,
                      const scs_float sc);
scs_float calc_norm_diff(const scs_float *a, const scs_float *b, scs_int l);
scs_float calc_norm_inf_diff(const scs_float *a, const scs_float *b, scs_int l);

#ifdef __cplusplus
}
#endif
#endif
