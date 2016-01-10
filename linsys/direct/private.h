#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include "glbopts.h"
#include "scs.h"
#include "cs.h"
#include "external/amd.h"
#include "external/ldl.h"
#include "../common.h"

struct PRIVATE_DATA {
    cs *L;         /* KKT, and factorization matrix L resp. */
    scs_float *D;  /* diagonal matrix of factorization */
    scs_int *P;    /* permutation of KKT matrix for factorization */
    scs_float *bp; /* workspace memory for solves */
    /* reporting */
    scs_float totalSolveTime;
};

#endif
