#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "linsys.h"
#include "scs_blas.h"

struct SCS_LIN_SYS_WORK {
  scs_int n;        /* number of QP variables */
  scs_int m;        /* number of QP constraints */
  scs_int n_plus_m; /* dimension of the linear system */

  /* Dense KKT matrix stored column-major */
  scs_float *kkt;     /* n_plus_m x n_plus_m dense matrix */
  blas_int *ipiv;     /* pivot indices from dgetrf */

  /* These are required for matrix updates */
  const ScsMatrix *A; /* does *not* own this memory */
  const ScsMatrix *P; /* does *not* own this memory */
  scs_float *diag_p;  /* Diagonal values of P */
};

#ifdef __cplusplus
}
#endif

#endif
