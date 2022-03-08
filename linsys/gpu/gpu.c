#include "gpu.h"

/* trans = CUSPARSE_OPERATION_TRANSPOSE is fastest for CSC */
void SCS(accum_by_a_gpu)(const ScsGpuMatrix *Ag,
                              cusparseOperation_t trans,
                              const cusparseDnVecDescr_t x,
                              cusparseDnVecDescr_t y,
                              cusparseHandle_t cusparse_handle,
                              size_t *buffer_size, void **buffer) {
  /* y += A*x or y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  size_t new_buffer_size = 0;
  CUDA_CHECK_ERR;

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, trans, &onef, Ag->descr, x,
   &onef, y, SCS_CUDA_FLOAT, SCS_MV_ALG, &new_buffer_size);

  if (new_buffer_size > *buffer_size) {
    if (*buffer != SCS_NULL) {
      cudaFree(*buffer);
    }
    cudaMalloc(buffer, *buffer_size);
    *buffer_size = new_buffer_size;
  }
  CUDA_CHECK_ERR;
  cusparseStatus_t status = CUSPARSE_GEN(SpMV)
  (cusparse_handle, trans, &onef, Ag->descr, x,
   &onef, y, SCS_CUDA_FLOAT, SCS_MV_ALG, buffer);
  CUDA_CHECK_ERR;
  scs_printf("status: %i\n", status);
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
  /* P is symmetric and FULL, so can use transpose here for speed */
  SCS(accum_by_a_gpu)(Pg, CUSPARSE_OPERATION_TRANSPOSE, x, y, cusparse_handle, buffer_size, buffer);
}

void SCS(free_gpu_matrix)(ScsGpuMatrix *A) {
  cudaFree(A->x);
  cudaFree(A->i);
  cudaFree(A->p);
  cusparseDestroySpMat(A->descr);
}
