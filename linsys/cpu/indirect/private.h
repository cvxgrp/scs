#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "csparse.h"
#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "scs_matrix.h"
#include "util.h" /* timer */
#include <string.h>

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
  /* eigCG deflation (Stathopoulos & Orginos 2010). During the deep
   * solve for g = K^{-1}h (cold, run to CG_BEST_TOL, once per metric
   * change) the CG recurrence is observed -- never perturbed -- and
   * approximate eigenvectors of the smallest eigenvalues of the
   * preconditioned operator are extracted by thick-restarted Lanczos on
   * CG's own scalars. Later warm solves Galerkin-project their initial
   * residual against this near-invariant subspace, which removes the
   * slowly-converging components; because the subspace is
   * near-invariant, CG does not reintroduce them. Discarded whenever
   * the metric changes; the deep re-solve re-harvests. */
  scs_float *dfl_p;  /* W: harvested eigenvectors, n x dfl_max */
  scs_float *dfl_Ap; /* A W images, n x dfl_max */
  scs_float *dfl_G;  /* Cholesky factor of W'AW, dfl_max x dfl_max */
  scs_float *dfl_c;  /* projection coefficients, dfl_max */
  scs_int dfl_max, dfl_count, dfl_harvest;
  /* eigCG window state (live only while harvesting) */
  scs_float *eig_V;    /* n x eig_win basis (Lanczos / locked Ritz) */
  scs_float *eig_T;    /* eig_win x eig_win projected operator */
  scs_float *eig_E;    /* eig_win x eig_win syev scratch */
  scs_float *eig_w;    /* eigenvalue scratch */
  scs_float *eig_S;    /* eig_win x 2 dfl_max restart stack */
  scs_float *eig_TS;   /* eig_win x 2 dfl_max scratch */
  scs_float *eig_H;    /* 2 dfl_max x 2 dfl_max small projection */
  scs_float *eig_G;    /* eig_win x 2 dfl_max restart transform */
  scs_float *eig_VS;   /* n x 2 dfl_max basis-update scratch */
  scs_float *eig_work; /* syev workspace */
  scs_int eig_win, eig_j, eig_jd, eig_step, eig_have_prev, eig_dead, eig_lwork;
  scs_float eig_pa, eig_pb; /* previous iteration's alpha / beta */
  /* reporting */
  scs_int tot_cg_its;
  const scs_float *diag_r;
};

#ifdef __cplusplus
}
#endif
#endif
