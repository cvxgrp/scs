#ifndef ACCEL_H_GUARD
#define ACCEL_H_GUARD

#include "glbopts.h"

typedef struct SCS_ACCEL Accel;

Accel* initAccel(Work * w);
void freeAccel(Accel * a);
scs_int accelerate(Work *w, scs_int iter);
char *getAccelSummary(const Info *info, Accel *a);

#endif
