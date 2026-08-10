#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "external/amd/amd.h"
#include "external/qdldl/qdldl.h"
#include "glbopts.h"
#include "linsys.h"
#include "scs_matrix.h"

struct SCS_LIN_SYS_WORK {
  scs_int m, n;       /* linear system dimensions */
  ScsMatrix *kkt, *L; /* KKT, and factorization matrix L resp. */
  scs_float *Dinv;    /* inverse diagonal matrix of factorization */
  scs_int *perm;      /* permutation of KKT matrix for factorization */
  scs_float *bp;      /* workspace memory for solves */
  scs_int *diag_r_idxs;
  scs_int soc_nnz;     /* (research) SOC border-column entry count */
  scs_int *soc_idxs;   /* their positions in kkt->x */
  scs_int nkkt;        /* KKT dimension: n + m + 2 * (SOC block count) */
  scs_float *bext;     /* bordered-rhs workspace (SOC only) */
  scs_int factorizations;
  /* ldl factorization workspace */
  scs_float *D, *fwork;
  scs_int *etree, *iwork, *Lnz;
  QDLDL_bool *bwork;
  scs_float *diag_p;
};

#ifdef __cplusplus
}
#endif
#endif
