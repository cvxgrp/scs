#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include "glbopts.h"
#include "scs.h"
#include "cs.h"
#include <math.h>
#include "../common.h"
#include "linAlg.h"

struct PRIVATE_DATA {
    scs_float *p; /* cg iterate  */
    scs_float *r; /* cg residual */
    scs_float *Gp;
    scs_float *tmp;
    AMatrix *At;
    /* preconditioning */
    scs_float *z;
    scs_float *M;
    /* reporting */
    scs_int totCgIts;
    scs_float totalSolveTime;
};

#endif
