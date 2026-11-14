#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <string.h>

/* A QP whose P is not positive semidefinite, P = diag(-10, 1), over the box
 * -1 <= x <= 1. The problem is not convex, and the direct linear solvers
 * that check the inertia of the KKT matrix at setup must refuse it, so scs
 * returns SCS_FAILED: QDLDL through its pivots, MKL Pardiso and cuDSS
 * through the number of positive eigenvalues, the dense solver through its
 * Cholesky factorization. What they test is whether R_x + P + A' R_y^{-1} A
 * is positive definite, so P11 = -10 is chosen to outweigh the A' R_y^{-1} A
 * term after equilibration by a wide margin; a mildly negative P11 such as
 * -0.1 would pass every backend's check. The indirect solvers and Apple
 * Accelerate do not check, so for them the test only runs a few iterations. */
static const char *test_non_psd_p(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *method = scs_get_lin_sys_method();
  int checks_inertia = strstr(method, "qdldl") || strstr(method, "mkl") ||
                       strstr(method, "cuDSS") || strstr(method, "cholesky");

  /* A = [I; -I], b = 1: x <= 1 and -x <= 1 */
  scs_float Ax[] = {1., -1., 1., -1.};
  scs_int Ai[] = {0, 2, 1, 3};
  scs_int Ap[] = {0, 2, 4};
  scs_float Px[] = {-10., 1.};
  scs_int Pi[] = {0, 1};
  scs_int Pp[] = {0, 1, 2};
  scs_float b[] = {1., 1., 1., 1.};
  scs_float c[] = {1., 1.};

  d->m = 4;
  d->n = 2;
  d->b = b;
  d->c = c;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = 4;
  d->A->n = 2;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;
  d->P->m = 2;
  d->P->n = 2;
  d->P->x = Px;
  d->P->i = Pi;
  d->P->p = Pp;
  k->l = 4;

  scs_set_default_settings(stgs);
  stgs->max_iters = 100;
  stgs->verbose = 0;

  exitflag = scs(d, k, stgs, sol, &info);
  scs_printf("non_psd_p: %s returned %li\n", method, (long)exitflag);
  if (checks_inertia) {
    mu_assert("non_psd_p: a direct solver accepted a P that is not PSD",
              exitflag == SCS_FAILED);
  }

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return 0;
}
