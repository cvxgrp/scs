#include "gpu.h"

void SCS(_accum_by_atrans_gpu)(const ScsGpuMatrix *Ag, const scs_float *x,
                               scs_float *y, cusparseHandle_t cusparse_handle,
                               size_t *buffer_size, void **buffer) {
  /* y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  cusparseDnVecDescr_t dn_vec_x = SCS_NULL, dn_vec_y = SCS_NULL;
  size_t new_buffer_size = 0;

  cusparseCreateDnVec(&dn_vec_x, Ag->m, (void *) x, SCS_CUDA_FLOAT);
  cusparseCreateDnVec(&dn_vec_y, Ag->n, (void *) y, SCS_CUDA_FLOAT);

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &onef, Ag->descr, dn_vec_x, &onef, dn_vec_y,
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
    &onef, Ag->descr, dn_vec_x, &onef, dn_vec_y,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    buffer);

  cusparseDestroyDnVec(dn_vec_x);
  cusparseDestroyDnVec(dn_vec_y);
}

void SCS(_accum_by_a_gpu)(const ScsGpuMatrix *Ag, const scs_float *x,
                          scs_float *y, cusparseHandle_t cusparse_handle,
                          size_t *buffer_size, void **buffer) {
  /* y += A*x
     x and y MUST be on GPU already
   */
  const scs_float onef = 1.0;
  cusparseDnVecDescr_t dn_vec_x = SCS_NULL, dn_vec_y = SCS_NULL;
  size_t new_buffer_size = 0;

  /* The A matrix idx pointers must be ORDERED */

  cusparseCreateDnVec(&dn_vec_x, Ag->n, (void *) x, SCS_CUDA_FLOAT);
  cusparseCreateDnVec(&dn_vec_y, Ag->m, (void *) y, SCS_CUDA_FLOAT);

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
    &onef, Ag->descr, dn_vec_x, &onef, dn_vec_y,
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
    &onef, Ag->descr, dn_vec_x, &onef, dn_vec_y,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    buffer);

  cusparseDestroyDnVec(dn_vec_x);
  cusparseDestroyDnVec(dn_vec_y);
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
