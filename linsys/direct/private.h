#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "external/amd/include/amd.h"
#include "external/qdldl/include/qdldl.h"
#include "glbopts.h"
#include "scs.h"

struct SCS_LIN_SYS_WORK {
  cs *L;            /* KKT, and factorization matrix L resp. */
  scs_float *Dinv;  /* inverse diagonal matrix of factorization */
  scs_int *P;       /* permutation of KKT matrix for factorization */
  scs_float *bp;    /* workspace memory for solves */
  scs_float total_solve_time; /* reporting */
};

#ifdef __cplusplus
}
#endif
#endif
