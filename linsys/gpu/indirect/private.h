#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "glbopts.h"
#include "gpu.h"
#include "linalg.h"
#include "scs.h"

struct SCS_LIN_SYS_WORK {
  scs_int n, m; /* linear system dimensions */
  /* reporting */
  scs_int tot_cg_its;
  scs_float *M; /* preconditioner on cpu */
  /* ALL BELOW HOSTED ON THE GPU */
  scs_float *p;       /* cg iterate, n  */
  scs_float *r;       /* cg residual, n */
  scs_float *Gp;      /* G * p, n */
  scs_float *bg;      /* b, n */
  scs_float *tmp_m;   /* m, used in mat_vec */
  scs_float *z;       /* preconditioned */
  scs_float *M_gpu;   /* preconditioner */
  const ScsMatrix *A; /* does *not* own this memory */
  const ScsMatrix *P; /* does *not* own this memory */
  ScsGpuMatrix *Ag;   /* A matrix on GPU */
  ScsGpuMatrix *Agt;  /* A trans matrix on GPU */
  ScsGpuMatrix *Pg;   /* P matrix on GPU */
  /* CUDA */
  cublasHandle_t cublas_handle;
  cusparseHandle_t cusparse_handle;
  /* CUSPARSE */
  size_t buffer_size;
  void *buffer;
  cusparseDnVecDescr_t dn_vec_m;   /* Dense vector of length m */
  cusparseDnVecDescr_t dn_vec_n;   /* Dense vector of length n */
  cusparseDnVecDescr_t dn_vec_n_p; /* Dense vector of length n */

  /* rho terms */
  scs_float *r_x_gpu;
  scs_float *inv_r_y;     /* inverse R_y */
  scs_float *inv_r_y_gpu; /* inverse R_y on GPU */
};

#ifdef __cplusplus
}
#endif
#endif
