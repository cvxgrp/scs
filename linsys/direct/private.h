#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "cs.h"
#include "external/amd.h"
#include "external/ldl.h"
#include "glbopts.h"
#include "scs.h"

struct SCS_LIN_SYS_WORK {
  cs *L;         /* KKT, and factorization matrix L resp. */
  scs_float *D;  /* diagonal matrix of factorization */
  scs_int *P;    /* permutation of KKT matrix for factorization */
  scs_float *bp; /* workspace memory for solves */
  /* reporting */
  scs_float total_solve_time;
};

#ifdef __cplusplus
}
#endif
#endif
