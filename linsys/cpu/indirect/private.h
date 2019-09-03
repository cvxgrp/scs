#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "amatrix.h"
#include "glbopts.h"
#include "linalg.h"
#include "scs.h"

struct SCS_LIN_SYS_WORK {
  scs_float *p; /* cg iterate  */
  scs_float *r; /* cg residual */
  scs_float *Gp;
  scs_float *tmp;
  ScsMatrix *At;
  /* preconditioning */
  scs_float *z;
  scs_float *M;
  /* reporting */
  scs_int tot_cg_its;
  scs_float total_solve_time;
};

#ifdef __cplusplus
}
#endif
#endif
