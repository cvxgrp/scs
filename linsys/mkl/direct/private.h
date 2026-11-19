#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "linsys.h"
#include "scs_matrix.h"

struct SCS_LIN_SYS_WORK {
  ScsMatrix *kkt; /* Upper triangular KKT matrix (in CSR format) */
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
  scs_int one_thread; /* Pardiso on one MKL thread (small system), see
                        PARDISO_MIN_THREADED_NNZ in private.c */

  /* These are required for matrix updates */
  scs_int *diag_r_idxs; /* indices where R appears */
  scs_float *diag_p;    /* Diagonal values of P */

  /* The multiplies by K of scs_solve_lin_sys (its warm start, refinement
   * and residual, see private.c), with K in MKL's sparse BLAS. prod_ok = 1
   * once set up; setup fails if they cannot be. prod_none = 1 when, by the
   * size rule in prod_setup or in a build with MKL_NO_OWN_PRODUCTS, this
   * solver does not compute them and they are SCS's products
   * (linsys/scs_matrix.h). */
  scs_int prod_ok;
  scs_int prod_none;
  void *prod_k;          /* MKL sparse_matrix_t handle of K */
  void *prod_ind;        /* column indices in MKL's integer type, or NULL
                            when kkt->i already has that type */
  void *prod_ptr;        /* row pointers likewise */

  /* the warm start, refinement and residual of scs_solve_lin_sys
   * (linsys/scs_matrix.h), with prod_residual unless prod_none */
  ScsLinSysShared sh;
};

#ifdef __cplusplus
}
#endif

#endif
