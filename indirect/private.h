#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cs.h"
#include "scs.h"
#include <math.h>

struct PRIVATE_DATA{
  double *p;  // cg iterate
  double *r;  // cg residual
  double *x;
  double * Ap;
  double * tmp;

  double * Atx;
  int * Ati;
  int * Atp;
};


#endif
