#include "util_spectral_cones.h"

// checks if all entries of a vector are nonnegative
bool is_non_neg(const scs_float *x, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        if (x[i] < 0.0)
        {
            return false;
        }
    }
    return true;
}

bool is_pos(const scs_float *x, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        if (x[i] <= 0.0)
        {
            return false;
        }
    }
    return true;
}

bool is_negative(const scs_float *x, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        if (x[i] >= 0.0)
        {
            return false;
        }
    }
    return true;
}

void non_neg_proj(const scs_float *src, scs_float *dst, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        dst[i] = (src[i] > 0.0) ? src[i] : 0.0;
    }
}

void print_vector(const scs_float *x, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        printf("%f ", x[i]);
    }
    printf("\n");
}

scs_float min_vec(const scs_float *vec, size_t n)
{
    scs_float minVal = vec[0];

    for (size_t i = 1; i < n; ++i)
    {
        if (vec[i] < minVal)
        {
            minVal = vec[i];
        }
    }

    return minVal;
}


scs_float sum_log(const scs_float *x, size_t n)
{
    scs_float sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sum += log(x[i]);
    }
    return sum;
}

