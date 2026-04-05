/*
 * Problem data normalization (equilibration) routines.
 * Applied before solving to improve numerical conditioning; undone
 * after solving to return the solution in the original scale.
 */

#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs_work.h"

void SCS(normalize_b_c)(ScsScaling *scal, scs_float *b, scs_float *c);
void SCS(normalize_sol)(ScsScaling *scal, ScsSolution *sol);
void SCS(un_normalize_sol)(ScsScaling *scal, ScsSolution *sol);

#ifdef __cplusplus
}
#endif
#endif
