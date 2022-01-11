#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "linalg.h"
#include "scs.h"
#include "scs_matrix.h"
#include <math.h>

struct SCS_LIN_SYS_WORK {
  scs_int n, m; /* linear system dimensions */
  scs_float *p; /* cg iterate  */
  scs_float *r; /* cg residual */
  scs_float *Gp;
  scs_float *tmp;
  const ScsMatrix *A; /* does *not* own this memory */
  const ScsMatrix *P; /* does *not* own this memory */
  ScsMatrix *At;      /* does own this memory */
  /* preconditioning */
  scs_float *z;
  scs_float *M;
  /* reporting */
  scs_int tot_cg_its;
  const scs_float *diag_r;
};

#ifdef __cplusplus
}
#endif
#endif
