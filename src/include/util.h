#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#include <stdlib.h>
#include <stdio.h>
#include "scs.h"

void tic(void); 
double toc(void); 
double tocq(void); 
void printConeData(Cone * k);
void printData(Data * d);
void printAll(Data * d, Work * w);


#endif
