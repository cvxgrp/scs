#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#include "glbopts.h"
#include "cones.h"
#include "amatrix.h"

void _accumByAtrans(scs_int n, scs_float * Ax, scs_int * Ai, scs_int * Ap, const scs_float *x, scs_float *y);
void _accumByA(scs_int n, scs_float * Ax, scs_int * Ai, scs_int * Ap, const scs_float *x, scs_float *y);

#endif
