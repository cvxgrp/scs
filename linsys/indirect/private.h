#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cs.h"
#include "scs.h"
#include <math.h>

struct PRIVATE_DATA{
  pfloat *p;  /* cg iterate  */
  pfloat *r;  /* cg residual */
  pfloat *x;
  pfloat * Ap;
  pfloat * tmp;

  pfloat * Atx;
  idxint * Ati;
  idxint * Atp;
};

#ifndef POWF
    #ifndef FLOAT
        #define POWF powf
    #else
        #define POWF pow
    #endif
#endif

#endif
