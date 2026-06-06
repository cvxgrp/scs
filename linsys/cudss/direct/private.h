#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "linsys.h"
#include <cuda_runtime.h>
#include <cudss.h>

/* cuDSS 0.8.0 renamed cudaDataType_t to cudssDataType_t (so CUDA_R_* enum
 * values become CUDSS_R_*) and added an offsetType parameter to
 * cudssMatrixCreateCsr. Detect via CUDSS_VER_MAJOR / CUDSS_VER_MINOR
 * (defined in <cudss.h>) and pick the right type tokens here so the call
 * sites stay readable. */
#if defined(CUDSS_VER_MAJOR) &&                                                \
    ((CUDSS_VER_MAJOR > 0) ||                                                  \
     (CUDSS_VER_MAJOR == 0 && CUDSS_VER_MINOR >= 8))
#define SCS_CUDSS_NEW_API 1
#else
#define SCS_CUDSS_NEW_API 0
#endif

#if SCS_CUDSS_NEW_API
#ifndef SFLOAT
#define SCS_CUDA_FLOAT CUDSS_R_64F
#else
#define SCS_CUDA_FLOAT CUDSS_R_32F
#endif
#ifndef DLONG
#define SCS_CUDA_INDEX CUDSS_R_32I
#else
#define SCS_CUDA_INDEX CUDSS_R_64I
#endif
#else
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
#endif

struct SCS_LIN_SYS_WORK {
  /* General problem dimensions */
  scs_int n;        /* number of QP variables */
  scs_int m;        /* number of QP constraints */
  scs_int n_plus_m; /* dimension of the linear system */

  /* CPU matrices and vectors */
  ScsMatrix *kkt; /* KKT matrix in CSR format */

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

  /* Pinned host memory for faster H<->D transfers */
  scs_float *h_b_pinned;   /* pinned host staging buffer for RHS */
  scs_float *h_sol_pinned; /* pinned host staging buffer for solution */

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
