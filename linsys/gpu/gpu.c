#include "gpu.h"

void SCS(accum_by_atrans_gpu)(const ScsGpuMatrix *Ag,
                              const cusparseDnVecDescr_t x,
                              cusparseDnVecDescr_t y,
                              cusparseHandle_t cusparse_handle,
                              size_t *buffer_size, void **buffer) {
  /* y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  size_t new_buffer_size = 0;

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &onef, Ag->descr, x,
   &onef, y, SCS_CUDA_FLOAT, SCS_CSRMV_ALG, &new_buffer_size);

  if (new_buffer_size > *buffer_size) {
    if (*buffer != SCS_NULL) {
      cudaFree(*buffer);
    }
    cudaMalloc(buffer, new_buffer_size);
    *buffer_size = new_buffer_size;
  }

  CUSPARSE_GEN(SpMV)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &onef, Ag->descr, x,
   &onef, y, SCS_CUDA_FLOAT, SCS_CSRMV_ALG, *buffer);
}

/* this is slow, use trans routine if possible */
void SCS(accum_by_a_gpu)(const ScsGpuMatrix *Ag, const cusparseDnVecDescr_t x,
                         cusparseDnVecDescr_t y,
                         cusparseHandle_t cusparse_handle, size_t *buffer_size,
                         void **buffer) {
  /* y += A*x
     x and y MUST be on GPU already
   */
  const scs_float onef = 1.0;
  size_t new_buffer_size = 0;

  /* The A matrix idx pointers must be ORDERED */
  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE, &onef, Ag->descr, x, &onef, y,
   SCS_CUDA_FLOAT, SCS_CSRMV_ALG, &new_buffer_size);

  if (new_buffer_size > *buffer_size) {
    if (*buffer != SCS_NULL) {
      cudaFree(*buffer);
    }
    cudaMalloc(buffer, new_buffer_size);
    *buffer_size = new_buffer_size;
  }

  CUSPARSE_GEN(SpMV)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE, &onef, Ag->descr, x, &onef, y,
   SCS_CUDA_FLOAT, SCS_CSRMV_ALG, *buffer);
}

/* This assumes that P has been made full (ie not triangular) and uses the
 * fact that the GPU is faster for general sparse matrices than for symmetric
 */
/* y += P*x
   x and y MUST be on GPU already
 */
void SCS(accum_by_p_gpu)(const ScsGpuMatrix *Pg, const cusparseDnVecDescr_t x,
                         cusparseDnVecDescr_t y,
                         cusparseHandle_t cusparse_handle, size_t *buffer_size,
                         void **buffer) {
  SCS(accum_by_atrans_gpu)(Pg, x, y, cusparse_handle, buffer_size, buffer);
}

void SCS(free_gpu_matrix)(ScsGpuMatrix *A) {
  cudaFree(A->x);
  cudaFree(A->i);
  cudaFree(A->p);
  cusparseDestroySpMat(A->descr);
}
