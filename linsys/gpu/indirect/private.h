#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "glbopts.h"
#include "linalg.h"
#include "scs.h"

#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

struct SCS_LIN_SYS_WORK {
  scs_int Annz; /* num non-zeros in A matrix */
  /* reporting */
  scs_int tot_cg_its;
  scs_float total_solve_time;
  /* CUDA */
  cublasHandle_t cublas_handle;
  cusparseHandle_t cusparse_handle;
  cusparseMatDescr_t descr;
  /* ALL BELOW HOSTED ON THE GPU */
  scs_float *p;     /* cg iterate, n  */
  scs_float *r;     /* cg residual, n */
  scs_float *Gp;    /* G * p, n */
  scs_float *bg;    /* b, n */
  scs_float *tmp_m; /* m, used in mat_vec */
  scs_float *z;     /* preconditioned */
  scs_float *M;     /* preconditioner */
  ScsMatrix *Ag;    /* A matrix on GPU */
  ScsMatrix *Agt;   /* A trans matrix on GPU */
};

#ifdef __cplusplus
}
#endif
#endif
