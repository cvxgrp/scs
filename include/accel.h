#ifndef ACCEL_H_GUARD
#define ACCEL_H_GUARD

#include "scs.h"
#include "glbopts.h"

ScsAccelWork *init_accel(ScsWork *w);
void free_accel(ScsAccelWork *a);
scs_int accelerate(ScsWork *w, scs_int iter);
char *get_accel_summary(const ScsInfo *info, ScsAccelWork *a);

#endif
