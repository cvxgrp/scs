#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C"
{
#endif

#include "csparse.h"
#include "linsys.h"
#include <hip/hip_runtime.h>
#include <hipsparse/hipsparse.h>

  struct SCS_LIN_SYS_WORK
  {
    // Host:
    ScsMatrix *kkt; /* Upper triangular KKT matrix (in CSR format) */
    scs_int n;      /* number of QP variables */
    scs_int m;      /* number of QP constraints */

    hipsparseHandle_t handle;  // HIPSPARSE handle
    hipsparseMatDescr_t descr; // Matrix descriptor

    // kkt matrix data on the device
    scs_float *d_vals;   // Non-zero values of the matrix (on device)
    scs_int *d_row_ptrs; // Row pointers (on device)
    scs_int *d_col_inds; // Column indices (on device)

    // Vectors for solving system Ax = b
    scs_float *d_b; // RHS vector b (on device)
    scs_float *d_x; // Solution vector x (on device)

    // LU decomposition info for the lower triangular matrix
    csrsv2Info_t info_LU; // Lower triangular solve info

    // Buffer for LU decomposition and solving
    void *buffer;       // Buffer for LU factorization and solving
    scs_int bufferSize; // Size of the buffer

    /* These are required for matrix updates */
    scs_int *diag_r_idxs; /* indices where R appears */
    scs_float *diag_p;    /* Diagonal values of P */
  };

#ifdef __cplusplus
}
#endif

#endif
