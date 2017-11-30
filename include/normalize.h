#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

void SCS(normalize_b_c)(ScsWork *w);
void SCS(calc_scaled_resids)(ScsWork *w, ScsResiduals *r);
void SCS(normalize_warm_start)(ScsWork *w);
void SCS(un_normalize_sol)(ScsWork *w, ScsSolution *sol);

#ifdef __cplusplus
}
#endif
#endif
