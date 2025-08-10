#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#ifndef SFLOAT
#define SCS_CUDA_FLOAT CUDA_R_64F
#else
#define SCS_CUDA_FLOAT CUDA_R_32F
#endif

#ifndef DLONG
#define SCS_CUDA_INDEX CUDA_R_32I
#else
#define SCS_CUDA_INDEX CUDA_R_64I
#endif

#include "csparse.h"
#include "linsys.h"
#include <cuda_runtime.h>
#include <cudss.h>

struct SCS_LIN_SYS_WORK {
  /* General problem dimensions */
  scs_int n;        /* number of QP variables */
  scs_int m;        /* number of QP constraints */
  scs_int n_plus_m; /* dimension of the linear system */

  /* CPU matrices and vectors */
  ScsMatrix *kkt; /* KKT matrix in CSR format */
  scs_float *sol; /* solution to the KKT system */

  /* cuDSS handle and descriptors */
  cudssHandle_t handle;    /* cuDSS library handle */
  cudssMatrix_t d_kkt_mat; /* cuDSS matrix descriptors */
  cudssMatrix_t d_b_mat;
  cudssMatrix_t d_sol_mat;

  /* Device memory for KKT matrix */
  scs_float *d_kkt_val;   /* device copy of KKT values */
  scs_int *d_kkt_row_ptr; /* device copy of KKT row pointers */
  scs_int *d_kkt_col_ind; /* device copy of KKT column indices */

  /* Device memory for vectors */
  scs_float *d_b;   /* device copy of right-hand side */
  scs_float *d_sol; /* device copy of solution */

  /* These are required for matrix updates */
  scs_int *diag_r_idxs; /* indices where R appears in the KKT matrix */
  scs_float *diag_p;    /* Diagonal values of P */

  /* cuDSS configuration */
  cudssConfig_t solver_config; /* cuDSS solver handle */
  cudssData_t solver_data;     /* cuDSS data handle */
};

#ifdef __cplusplus
}
#endif

#endif
