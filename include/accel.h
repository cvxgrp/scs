#ifndef ACCEL_H_GUARD
#define ACCEL_H_GUARD

#include "scs.h"
#include "glbopts.h"

ScsAccelWork *SCS(init_accel)(ScsWork *w);
void SCS(free_accel)(ScsAccelWork *a);
scs_int SCS(accelerate)(ScsWork *w, scs_int iter);
char *SCS(get_accel_summary)(const ScsInfo *info, ScsAccelWork *a);

#endif
