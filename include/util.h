#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#include <stdlib.h>
#include <stdio.h>
#include "scs.h"
#include "cones.h"

void tic(void);
pfloat toc(void);
pfloat tocq(void);
void printConeData(Cone * k);
void printData(Data * d);
void printWork(Data * d, Work * w);

#endif
