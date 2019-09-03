#include "gpu.h"

void SCS(_accum_by_atrans_gpu)(const ScsGpuMatrix *Ag, const scs_float *x,
                               scs_float *y, cusparseHandle_t cusparse_handle) {
  /* y += A'*x
     x and y MUST be on GPU already
  */
  const scs_float onef = 1.0;
  CUSPARSE(csrmv)
  (cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Ag->n, Ag->m, Ag->Annz,
   &onef, Ag->descr, Ag->x, Ag->p, Ag->i, x, &onef, y);
}

void SCS(_accum_by_a_gpu)(const ScsGpuMatrix *Ag, const scs_float *x,
                          scs_float *y, cusparseHandle_t cusparse_handle) {
  /* y += A*x
     x and y MUST be on GPU already
   */
  const scs_float onef = 1.0;
  /* The A matrix idx pointers must be ORDERED */
  CUSPARSE(csrmv)
  (cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE, Ag->n, Ag->m, Ag->Annz, &onef,
   Ag->descr, Ag->x, Ag->p, Ag->i, x, &onef, y);
}

void SCS(free_gpu_matrix)(ScsGpuMatrix *A) {
  cudaFree(A->x);
  cudaFree(A->i);
  cudaFree(A->p);
  cusparseDestroyMatDescr(A->descr);
}

void SCS(normalize_a)(ScsMatrix *A, const ScsSettings *stgs, const ScsCone *k,
                      ScsScaling *scal) {
  SCS(_normalize_a)(A, stgs, k, scal);
}

void SCS(un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                         const ScsScaling *scal) {
  SCS(_un_normalize_a)(A, stgs, scal);
}
