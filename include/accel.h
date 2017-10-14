#ifndef ACCEL_H_GUARD
#define ACCEL_H_GUARD

#include "glbopts.h"

typedef struct SCS_ACCEL Accel;

Accel* init_accel(Work * w);
void free_accel(Accel * a);
scs_int accelerate(Work *w, scs_int iter);
char *get_accel_summary(const Info *info, Accel *a);

#endif
