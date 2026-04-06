#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "linsys.h"
#include "scs_blas.h"

struct SCS_LIN_SYS_WORK {
  scs_int n; /* number of QP variables */
  scs_int m; /* number of QP constraints */

  scs_float *A_dense; /* dense A matrix, m x n column-major */
  scs_float *G;       /* Gram matrix R_x + P + A' R_y^{-1} A, n x n col-major */
  scs_float *r_y_inv; /* 1 / R_y diagonal, length m */
  scs_float *tmp_m;   /* workspace of length m */

  /* These are required for matrix updates */
  const ScsMatrix *P; /* does *not* own this memory */
  scs_float *diag_p;  /* Diagonal values of P */
};

#ifdef __cplusplus
}
#endif

#endif
