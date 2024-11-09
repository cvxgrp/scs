#ifndef UTILSPECTRALCONES_H
#define UTILSPECTRALCONES_H
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "scs_types.h"
#include "scs_blas.h"

#define IN_CONE -1
#define IN_NEGATIVE_DUAL_CONE -2
#define ANALYTICAL_SOL -3

// scs_float BLAS(dot)(const blas_int *n, const scs_float *x, const blas_int *incx,
//                     const scs_float *y, const blas_int *incy);

// void BLAS(axpy)(blas_int *n, const scs_float *a, const scs_float *x,
//                 blas_int *incx, scs_float *y, blas_int *incy);


long get_time_in_microseconds();
bool is_non_neg(const scs_float *x, size_t n);
bool is_pos(const scs_float *x, size_t n);
bool is_negative(const scs_float *x, size_t n);
void non_neg_proj(const scs_float *src, scs_float *dst, size_t n);
void print_vector(const scs_float *x, size_t n);
scs_float min_vec(const scs_float *vec, size_t n);
scs_float sum_log(const scs_float *x, size_t n);


// used for sorting in ell1-norm cone and sum of largest cone.
typedef struct
{
    scs_float value;
    int index;
} Value_index;

typedef struct 
{
    int iter;
    
    // if plain Newton computed the projection
    int NewtonSuccess; 
    
    // dual_res, pri_res, complementarity for the projection problem
    scs_float residuals[3];
} newton_stats;

#endif