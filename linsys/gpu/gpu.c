#include "gpu.h"

void SCS(_accum_by_atrans_gpu)(const ScsGpuMatrix *Ag, const scs_float *x,
                               scs_float *y, cusparseHandle_t cusparse_handle) {
  /* y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  cusparseDnVecDescr_t dnVecX = SCS_NULL, dnVecY = SCS_NULL;
  size_t bufferSize = 0;
  void *tmpBuffer = SCS_NULL;

  cusparseCreateDnVec(&dnVecX, Ag->m, (void *) x, SCS_CUDA_FLOAT);
  cusparseCreateDnVec(&dnVecY, Ag->n, (void *) y, SCS_CUDA_FLOAT);

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &onef, Ag->descr, dnVecX, &onef, dnVecY,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    &bufferSize);
  cudaMalloc(&tmpBuffer, bufferSize); 
  CUSPARSE_GEN(SpMV)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    &onef, Ag->descr, dnVecX, &onef, dnVecY,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    tmpBuffer);

  cusparseDestroyDnVec(dnVecX);
  cusparseDestroyDnVec(dnVecY);
  cudaFree(tmpBuffer);
}

void SCS(_accum_by_a_gpu)(const ScsGpuMatrix *Ag, const scs_float *x,
                          scs_float *y, cusparseHandle_t cusparse_handle) {
  /* y += A*x
     x and y MUST be on GPU already
   */
  const scs_float onef = 1.0;
  cusparseDnVecDescr_t dnVecX = SCS_NULL, dnVecY = SCS_NULL;
  size_t bufferSize = 0;
  void *tmpBuffer = SCS_NULL;

  /* The A matrix idx pointers must be ORDERED */

  cusparseCreateDnVec(&dnVecX, Ag->n, (void *) x, SCS_CUDA_FLOAT);
  cusparseCreateDnVec(&dnVecY, Ag->m, (void *) y, SCS_CUDA_FLOAT);

  CUSPARSE_GEN(SpMV_bufferSize)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
    &onef, Ag->descr, dnVecX, &onef, dnVecY,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    &bufferSize);
  cudaMalloc(&tmpBuffer, bufferSize);
  CUSPARSE_GEN(SpMV)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
    &onef, Ag->descr, dnVecX, &onef, dnVecY,
    SCS_CUDA_FLOAT, SCS_CSRMV_ALG,
    tmpBuffer);

  cusparseDestroyDnVec(dnVecX);
  cusparseDestroyDnVec(dnVecY);
  cudaFree(tmpBuffer);
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
