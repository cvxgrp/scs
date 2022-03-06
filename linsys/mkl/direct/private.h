#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "linsys.h"

struct SCS_LIN_SYS_WORK {
  ScsMatrix *kkt; /* Upper triangular KKT matrix (in CSR format) */
  scs_float *sol; /* solution to the KKT system */
  scs_int n;      /* number of QP variables */
  scs_int m;      /* number of QP constraints */

  /* Pardiso variables */
  void *pt[64];      /* internal solver memory pointer pt */
  scs_int iparm[64]; /* Pardiso control parameters */
  scs_int n_plus_m;  /* dimension of the linear system */
  scs_int mtype;     /* matrix type (-2 for real and symmetric indefinite) */
  scs_int nrhs;      /* number of right-hand sides (1) */
  scs_int maxfct;    /* maximum number of factors (1) */
  scs_int mnum;      /* indicates matrix for the solution phase (1) */
  scs_int phase;     /* control the execution phases of the solver */
  scs_int error;     /* the error indicator (0 for no error) */
  scs_int msglvl;    /* Message level information (0 for no output) */

  /* These are required for matrix updates */
  scs_int *diag_r_idxs; /* indices where R appears */
  scs_float *diag_p;    /* Diagonal values of P */
};

#ifdef __cplusplus
}
#endif

#endif
