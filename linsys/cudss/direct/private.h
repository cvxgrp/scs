#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "linsys.h"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cudss.h>
#include <cusparse.h>

/* cuDSS 0.8.0 renamed cudaDataType_t to cudssDataType_t (so CUDA_R_* enum
 * values become CUDSS_R_*) and added an offsetType parameter to
 * cudssMatrixCreateCsr. cudss.h exposes CUDSS_VERSION = MAJOR*10000 +
 * MINOR*100 + PATCH; gate on that so the same source still builds against
 * older cuDSS (<=0.7.x). */
#if defined(CUDSS_VERSION) && CUDSS_VERSION >= 800
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

/* value and index types of the cuSPARSE product (scs_lin_sys_residual) */
#ifndef SFLOAT
#define SCS_CUSPARSE_FLOAT CUDA_R_64F
#else
#define SCS_CUSPARSE_FLOAT CUDA_R_32F
#endif
#ifndef DLONG
#define SCS_CUSPARSE_INDEX CUSPARSE_INDEX_32I
#else
#define SCS_CUSPARSE_INDEX CUSPARSE_INDEX_64I
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

  /* Products with K for SCS's per-iteration refinement (scs_lin_sys_residual,
   * scs_lin_sys_keep, scs_lin_sys_correct), on the GPU. cuSPARSE has no
   * symmetric SpMV, so K is stored a second time with both triangles (CSR,
   * columns sorted in each row); its values are kkt->x[full_src[k]] and
   * follow scs_update_lin_sys_diag_r. The guess x0 and the residual r stay
   * on the GPU. products_ok = 1 when K has at least CUDSS_MIN_PRODUCT_NNZ
   * nonzeros (setup fails if the products cannot be set up then), else 0:
   * scs_lin_sys_residual returns "not provided" and SCS computes them. */
  scs_int products_ok;
  scs_int full_nnz;        /* nonzeros of K with both triangles */
  scs_int *full_src;       /* host: index in kkt->x of each value */
  scs_float *full_val_h;   /* host staging of the values */
  scs_float *d_full_val;   /* device CSR of K, both triangles */
  scs_int *d_full_row_ptr;
  scs_int *d_full_col_ind;
  scs_float *d_r;          /* device r = b - K x */
  scs_float *d_x0;         /* device guess x0 (scs_lin_sys_keep) */
  scs_float *d_xv;         /* device copy of an x given to the residual */
  void *d_spmv_buf;        /* cuSPARSE SpMV workspace */
  cusparseHandle_t sp_handle;
  cublasHandle_t bl_handle;
  cusparseSpMatDescr_t sp_k;
  cusparseDnVecDescr_t sp_r, sp_x0, sp_xv;
};

#ifdef __cplusplus
}
#endif

#endif
