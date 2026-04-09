#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include <Accelerate/Accelerate.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "glbopts.h"
#include "linsys.h"
#include "scs_matrix.h"

struct SCS_LIN_SYS_WORK {
  scs_int m, n;
  ScsMatrix *kkt;         /* Upper triangular KKT matrix (CSC format) */
  scs_int *diag_r_idxs;  /* Indices where R appears in kkt->x */
  scs_float *diag_p;     /* Diagonal values of P */
  long *col_starts;       /* Column starts converted to long for Accelerate */
#ifdef SFLOAT
  SparseMatrix_Float accel_mat;
  SparseOpaqueFactorization_Float factorization;
#else
  SparseMatrix_Double accel_mat;
  SparseOpaqueFactorization_Double factorization;
#endif
  scs_int factorizations; /* Number of successful factorizations */
};

#ifdef __cplusplus
}
#endif
#endif
