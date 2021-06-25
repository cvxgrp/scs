#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

void SCS(normalize_sol)(ScsWork *w, ScsSolution *sol);
void SCS(un_normalize_sol)(ScsWork *w, ScsSolution *sol);
void SCS(un_normalize_primal)(ScsWork *w, scs_float *r);
void SCS(un_normalize_dual)(ScsWork *w, scs_float *r);

#ifdef __cplusplus
}
#endif
#endif
