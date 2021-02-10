#include "gpu.h"

void SCS(_accum_by_atrans_gpu)(const ScsGpuMatrix *Ag, const cusparseDnVecDescr_t x,
                               cusparseDnVecDescr_t y, cusparseHandle_t cusparse_handle,
                               size_t *buffer_size, void **buffer) {
  /* y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  size_t new_buffer_size = 0;

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &onef, Ag->descr, x, &onef, y,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    &new_buffer_size);

  if (new_buffer_size > *buffer_size) {
    if (*buffer != SCS_NULL) {
      cudaFree(*buffer);
    }
    cudaMalloc(buffer, *buffer_size);
    *buffer_size = new_buffer_size;
  }

  CUSPARSE_GEN(SpMV)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &onef, Ag->descr, x, &onef, y,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    buffer);
}

void SCS(_accum_by_a_gpu)(const ScsGpuMatrix *Ag, const cusparseDnVecDescr_t x,
                          cusparseDnVecDescr_t y, cusparseHandle_t cusparse_handle,
                          size_t *buffer_size, void **buffer) {
  /* y += A*x
     x and y MUST be on GPU already
   */
  const scs_float onef = 1.0;
  size_t new_buffer_size = 0;

  /* The A matrix idx pointers must be ORDERED */

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
    &onef, Ag->descr, x, &onef, y,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    &new_buffer_size);

  if (new_buffer_size > *buffer_size) {
    if (*buffer != SCS_NULL) {
      cudaFree(*buffer);
    }
    cudaMalloc(buffer, *buffer_size);
    *buffer_size = new_buffer_size;
  }

  CUSPARSE_GEN(SpMV)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
    &onef, Ag->descr, x, &onef, y,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    buffer);
}

void SCS(free_gpu_matrix)(ScsGpuMatrix *A) {
  cudaFree(A->x);
  cudaFree(A->i);
  cudaFree(A->p);
  cusparseDestroySpMat(A->descr);
}

void SCS(normalize_a)(ScsMatrix *A, const ScsSettings *stgs, const ScsCone *k,
                      ScsScaling *scal) {
  SCS(_normalize_a)(A, stgs, k, scal);
}

void SCS(un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                         const ScsScaling *scal) {
  SCS(_un_normalize_a)(A, stgs, scal);
}
