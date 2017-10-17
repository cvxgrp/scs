#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#include "scs.h"

void normalize_b_c(ScsWork *w);
void calc_scaled_resids(ScsWork *w, ScsResiduals *r);
void normalize_warm_start(ScsWork *w);
void un_normalize_sol(ScsWork *w, ScsSolution *sol);

#endif
