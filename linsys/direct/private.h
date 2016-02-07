#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include "glbopts.h"
#include "scs.h"
#include "../common.h"

struct PRIVATE_DATA {
    scs_float *L;
    /* reporting */
    scs_float totalSolveTime;
};

#endif
