#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "amatrix.h"
#include "external/amd/amd.h"
#include "external/qdldl/qdldl.h"
#include "glbopts.h"
#include "scs.h"

typedef struct SPARSE_MATRIX _cs;
struct SCS_LIN_SYS_WORK {
  _cs *kkt, *L;            /* KKT, and factorization matrix L resp. */
  scs_float *Dinv;  /* inverse diagonal matrix of factorization */
  scs_int *perm;       /* permutation of KKT matrix for factorization */
  scs_float *bp;    /* workspace memory for solves */
  scs_int *scale_idxs;
};

#ifdef __cplusplus
}
#endif
#endif
