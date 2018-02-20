#include "linalg.h"
#include <math.h>

/* x = b*a */
void SCS(set_as_scaled_array)(scs_float *x, const scs_float *a,
                              const scs_float b, scs_int len) {
  scs_int i;
  for (i = 0; i < len; ++i) x[i] = b * a[i];
}

/* a *= b */
void SCS(scale_array)(scs_float *a, const scs_float b, scs_int len) {
  scs_int i;
  for (i = 0; i < len; ++i) a[i] *= b;
}

/* x'*y */
scs_float SCS(dot)(const scs_float *x, const scs_float *y, scs_int len) {
  scs_int i;
  scs_float ip = 0.0;
  for (i = 0; i < len; ++i) {
    ip += x[i] * y[i];
  }
  return ip;
}

/* ||v||_2^2 */
scs_float SCS(norm_sq)(const scs_float *v, scs_int len) {
  scs_int i;
  scs_float nmsq = 0.0;
  for (i = 0; i < len; ++i) {
    nmsq += v[i] * v[i];
  }
  return nmsq;
}

/* ||v||_2 */
scs_float SCS(norm)(const scs_float *v, scs_int len) {
  return SQRTF(SCS(norm_sq)(v, len));
}

scs_float SCS(norm_inf)(const scs_float *a, scs_int l) {
  scs_float tmp, max = 0.0;
  scs_int i;
  for (i = 0; i < l; ++i) {
    tmp = ABS(a[i]);
    if (tmp > max) {
      max = tmp;
    }
  }
  return max;
}

/* saxpy a += sc*b */
void SCS(add_scaled_array)(scs_float *a, const scs_float *b, scs_int n,
                           const scs_float sc) {
  scs_int i;
  for (i = 0; i < n; ++i) {
    a[i] += sc * b[i];
  }
}

scs_float SCS(norm_diff)(const scs_float *a, const scs_float *b, scs_int l) {
  scs_float nm_diff = 0.0, tmp;
  scs_int i;
  for (i = 0; i < l; ++i) {
    tmp = (a[i] - b[i]);
    nm_diff += tmp * tmp;
  }
  return SQRTF(nm_diff);
}

scs_float SCS(norm_inf_diff)(const scs_float *a, const scs_float *b,
                             scs_int l) {
  scs_float tmp, max = 0.0;
  scs_int i;
  for (i = 0; i < l; ++i) {
    tmp = ABS(a[i] - b[i]);
    if (tmp > max) {
      max = tmp;
    }
  }
  return max;
}
