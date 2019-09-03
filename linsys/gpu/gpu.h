#ifndef SCSGPU_H_GUARD
#define SCSGPU_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

#include "amatrix.h"
#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "scs.h"
#include "util.h"

#define CUDA_CHECK_ERR                                                    \
  do {                                                                    \
    cudaError_t err = cudaGetLastError();                                 \
    if (err != cudaSuccess) {                                             \
      printf("%s:%d:%s\n ERROR_CUDA: %s\n", __FILE__, __LINE__, __func__, \
             cudaGetErrorString(err));                                    \
    }                                                                     \
  } while (0)

#ifndef EXTRA_VERBOSE
#ifndef SFLOAT
#define CUBLAS(x) cublasD##x
#define CUSPARSE(x) cusparseD##x
#else
#define CUBLAS(x) cublasS##x
#define CUSPARSE(x) cusparseS##x
#endif
#else
#ifndef SFLOAT
#define CUBLAS(x) \
  CUDA_CHECK_ERR; \
  cublasD##x
#define CUSPARSE(x) \
  CUDA_CHECK_ERR;   \
  cusparseD##x
#else
#define CUBLAS(x) \
  CUDA_CHECK_ERR; \
  cublasS##x
#define CUSPARSE(x) \
  CUDA_CHECK_ERR;   \
  cusparseS##x
#endif
#endif

/*
 CUDA matrix routines only for CSR, not CSC matrices:
    CSC             CSR             GPU     Mult
    A (m x n)       A' (n x m)      Ag      accum_by_a_trans_gpu
    A'(n x m)       A  (m x n)      Agt     accum_by_a_gpu
*/

/* this struct defines the data matrix A on GPU */
typedef struct SCS_GPU_A_DATA_MATRIX {
  /* A is supplied in column compressed format */
  scs_float *x; /* A values, size: NNZ A */
  scs_int *i;   /* A row index, size: NNZ A */
  scs_int *p;   /* A column pointer, size: n+1 */
  scs_int m, n; /* m rows, n cols */
  scs_int Annz; /* num non-zeros in A matrix */
  /* CUDA */
  cusparseMatDescr_t descr;
} ScsGpuMatrix;

void SCS(_accum_by_atrans_gpu)(const ScsGpuMatrix *A, const scs_float *x,
                               scs_float *y, cusparseHandle_t cusparse_handle);

void SCS(_accum_by_a_gpu)(const ScsGpuMatrix *A, const scs_float *x,
                          scs_float *y, cusparseHandle_t cusparse_handle);

void SCS(free_gpu_matrix)(ScsGpuMatrix *A);

#ifdef __cplusplus
}
#endif
#endif
