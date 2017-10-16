#ifndef ACCEL_H_GUARD
#define ACCEL_H_GUARD

#include "glbopts.h"

typedef struct SCS_ACCEL_WORK ScsAccelWork;

ScsAccelWork* init_accel(ScsWork * w);
void free_accel(ScsAccelWork * a);
scs_int accelerate(ScsWork *w, scs_int iter);
char *get_accel_summary(const ScsInfo *info, ScsAccelWork *a);

#endif
