#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#include "glbopts.h"
#include "cones.h"
#include "amatrix.h"

/* contains routines common to direct and indirect sparse solvers */
scs_int validateLinSys(Data *d);
void normalizeA(Data * d, Work * w, Cone * k);
void unNormalizeA(Data *d, Work * w);
#endif
