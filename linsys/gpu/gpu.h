#ifndef SCS_GPU_H_GUARD
#define SCS_GPU_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: Do we need this?

#include <cuda.h>

*/

#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

#define CUDA_CHECK_ERR                                                         \
  do {                                                                         \
    cudaDeviceSynchronize();                                                   \
    cudaError_t err = cudaGetLastError();                                      \
    if (err != cudaSuccess) {                                                  \
      scs_printf("%s:%d:%s\n ERROR_CUDA (#): %s\n", __FILE__, __LINE__,        \
                 __func__, cudaGetErrorString(err));                           \
    }                                                                          \
  } while (0)

#if VERBOSITY == 0
#ifndef SFLOAT
#define CUBLAS(x) cublasD##x
#define CUBLASI(x) cublasId##x
#else
#define CUBLAS(x) cublasS##x
#define CUBLASI(x) cublasIs##x
#endif
#define CUSPARSE_GEN(x) cusparse##x
#else
#ifndef SFLOAT
#define CUBLAS(x)                                                              \
  CUDA_CHECK_ERR;                                                              \
  cublasD##x
#define CUBLASI(x)                                                             \
  CUDA_CHECK_ERR;                                                              \
  cublasId##x
#else
#define CUBLAS(x)                                                              \
  CUDA_CHECK_ERR;                                                              \
  cublasS##x
#define CUBLASI(x)                                                             \
  CUDA_CHECK_ERR;                                                              \
  cublasIs##x
#endif
#define CUSPARSE_GEN(x)                                                        \
  CUDA_CHECK_ERR;                                                              \
  cusparse##x
#endif

#ifndef SFLOAT
#define SCS_CUDA_FLOAT CUDA_R_64F
#else
#define SCS_CUDA_FLOAT CUDA_R_32F
#endif

#ifndef DLONG
#define SCS_CUSPARSE_INDEX CUSPARSE_INDEX_32I
#else
#define SCS_CUSPARSE_INDEX CUSPARSE_INDEX_64I
#endif

#define SCS_CSRMV_ALG CUSPARSE_SPMV_CSR_ALG1
#define SCS_CSR2CSC_ALG CUSPARSE_CSR2CSC_ALG1

/*
 CUDA matrix routines only for CSR, not CSC matrices:
    CSC             CSR             GPU     Mult
    A (m x n)       A' (n x m)      Ag      accum_by_a_trans_gpu
    A'(n x m)       A  (m x n)      Agt     accum_by_a_gpu
*/

/* this struct defines the data matrix on GPU */
typedef struct SCS_GPU_DATA_MATRIX {
  /* A is supplied in column compressed format */
  scs_float *x; /* values, size: NNZ */
  scs_int *i;   /* row index, size: NNZ */
  scs_int *p;   /* column pointer, size: n+1 */
  scs_int m, n; /* m rows, n cols */
  scs_int nnz;  /* num non-zeros in matrix */
  /* CUDA */
  cusparseSpMatDescr_t descr;
} ScsGpuMatrix;

void SCS(accum_by_atrans_gpu)(const ScsGpuMatrix *A,
                              const cusparseDnVecDescr_t x,
                              cusparseDnVecDescr_t y,
                              cusparseHandle_t cusparse_handle,
                              size_t *buffer_size, void **buffer);

void SCS(accum_by_a_gpu)(const ScsGpuMatrix *A, const cusparseDnVecDescr_t x,
                         cusparseDnVecDescr_t y,
                         cusparseHandle_t cusparse_handle, size_t *buffer_size,
                         void **buffer);

void SCS(accum_by_p_gpu)(const ScsGpuMatrix *P, const cusparseDnVecDescr_t x,
                         cusparseDnVecDescr_t y,
                         cusparseHandle_t cusparse_handle, size_t *buffer_size,
                         void **buffer);

void SCS(free_gpu_matrix)(ScsGpuMatrix *A);

#ifdef __cplusplus
}
#endif
#endif
