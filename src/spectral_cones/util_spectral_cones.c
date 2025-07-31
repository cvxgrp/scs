#include "util_spectral_cones.h"

bool is_pos(const scs_float *x, scs_int n) {
  for (scs_int i = 0; i < n; ++i) {
    if (x[i] <= 0.0) {
      return false;
    }
  }
  return true;
}

bool is_negative(const scs_float *x, scs_int n) {
  for (scs_int i = 0; i < n; ++i) {
    if (x[i] >= 0.0) {
      return false;
    }
  }
  return true;
}

void non_neg_proj(const scs_float *src, scs_float *dst, scs_int n) {
  for (scs_int i = 0; i < n; ++i) {
    dst[i] = (src[i] > 0.0) ? src[i] : 0.0;
  }
}

void print_vector(const scs_float *x, scs_int n) {
  for (scs_int i = 0; i < n; ++i) {
    printf("%f ", x[i]);
  }
  printf("\n");
}

scs_float min_vec(const scs_float *vec, scs_int n) {
  scs_float minVal = vec[0];

  for (scs_int i = 1; i < n; ++i) {
    if (vec[i] < minVal) {
      minVal = vec[i];
    }
  }

  return minVal;
}

scs_float sum_log(const scs_float *x, scs_int n) {
  scs_float sum = 0.0;
  for (scs_int i = 0; i < n; ++i) {
    sum += log(x[i]);
  }
  return sum;
}
