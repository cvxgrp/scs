#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "external/amd/amd.h"
#include "external/qdldl/qdldl.h"
#include "glbopts.h"
#include "scs.h"
#include "scs_matrix.h"

struct SCS_LIN_SYS_WORK {
  scs_int m, n;    /* linear system dimensions */
  csc *kkt, *L;    /* KKT, and factorization matrix L resp. */
  scs_float *Dinv; /* inverse diagonal matrix of factorization */
  scs_int *perm;   /* permutation of KKT matrix for factorization */
  scs_float *bp;   /* workspace memory for solves */
  scs_int *diag_r_idxs;
  scs_int factorizations;
  /* ldl factorization workspace */
  scs_float *D, *fwork;
  scs_int *etree, *iwork, *Lnz, *bwork;
  scs_float *diag_p;
};

#ifdef __cplusplus
}
#endif
#endif
