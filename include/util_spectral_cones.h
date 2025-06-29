#ifndef UTILSPECTRALCONES_H
#define UTILSPECTRALCONES_H

#ifdef __cplusplus
extern "C" {
#endif

#include "scs_blas.h"
#include "scs_types.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#define IN_CONE -1
#define IN_NEGATIVE_DUAL_CONE -2
#define ANALYTICAL_SOL -3

bool is_pos(const scs_float *x, scs_int n);
bool is_negative(const scs_float *x, scs_int n);
void non_neg_proj(const scs_float *src, scs_float *dst, scs_int n);
scs_float sum_log(const scs_float *x, scs_int n);
scs_float min_vec(const scs_float *vec, scs_int n);

// used for sorting in ell1-norm cone and sum of largest cone.
typedef struct {
  scs_float value;
  int index;
} Value_index;

typedef struct {
  int iter;

  // if plain Newton computed the projection or if an IPM was used
  int newton_success;

  // dual_res, pri_res, complementarity for the projection problem
  scs_float residuals[3];
} Newton_stats;

#ifdef __cplusplus
}
#endif
#endif
