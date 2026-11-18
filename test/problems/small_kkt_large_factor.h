/*
 * A QP whose KKT matrix is small but whose factor is not: min 1/2 ||x||^2
 * subject to L x = b, with L the 5-point Laplacian of a k x k grid (k = 50,
 * n = m = 2500). The KKT matrix has about 1.7e4 nonzeros in its upper
 * triangle; a fill-reducing ordering gives a factor of about 1.3e5 (AMD).
 *
 * With MKL Pardiso this takes the branch of PARDISO_MIN_THREADED_NNZ (5e4,
 * linsys/mkl/direct/private.c) where the symbolic phase first runs on one
 * thread (small KKT matrix) and is then redone on MKL's threads (large
 * factor), after releasing the Pardiso handle. The test checks that the
 * solve is still right after that. The other backends simply solve it.
 * Tolerances are SCS's 1e-5 and a residual check at 1e-4 of ||b||: the test
 * is about the handle logic, not accuracy. Skipped in single precision, which
 * the MKL handle logic does not depend on.
 */
#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

static const char *small_kkt_large_factor(void) {
#ifdef SFLOAT
  return mu_skipped;
#else
  const scs_int k = 50, n = k * k;
  ScsCone *cone = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsMatrix *P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsInfo info = {0};
  scs_int j, r, c, nz = 0, exitflag;
  scs_float *res, res_max = 0., b_max = 0.;
  const char *fail = SCS_NULL;

  /* A = L, column by column (L is symmetric, rows in increasing order) */
  A->m = n;
  A->n = n;
  A->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  A->i = (scs_int *)scs_calloc(5 * n, sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(5 * n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    r = j / k;
    c = j % k;
    if (r > 0) {
      A->i[nz] = j - k;
      A->x[nz++] = -1.;
    }
    if (c > 0) {
      A->i[nz] = j - 1;
      A->x[nz++] = -1.;
    }
    A->i[nz] = j;
    A->x[nz++] = 4.;
    if (c < k - 1) {
      A->i[nz] = j + 1;
      A->x[nz++] = -1.;
    }
    if (r < k - 1) {
      A->i[nz] = j + k;
      A->x[nz++] = -1.;
    }
    A->p[j + 1] = nz;
  }
  /* P = I (upper triangle) */
  P->m = n;
  P->n = n;
  P->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  P->i = (scs_int *)scs_calloc(n, sizeof(scs_int));
  P->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    P->i[j] = j;
    P->x[j] = 1.;
    P->p[j + 1] = j + 1;
  }
  d->m = n;
  d->n = n;
  d->A = A;
  d->P = P;
  d->b = (scs_float *)scs_calloc(n, sizeof(scs_float));
  d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    /* any right-hand side, scaled by h^2 = 1/k^2 so that x and the objective
     * are O(1): with b = O(1) the objective is about 1e7 and rounding in the
     * duality gap reaches verify_solution_correct's 1e-7 tolerance */
    d->b[j] = (1. + 0.5 * sin((scs_float)j)) / (scs_float)(k * k);
    b_max = MAX(b_max, ABS(d->b[j]));
  }
  cone->z = n; /* L x = b; L is nonsingular, so x = L^{-1} b */

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-5;
  stgs->eps_rel = 1e-5;
  stgs->verbose = 0;

  exitflag = scs(d, cone, stgs, sol, &info);
  mu_assert("small_kkt_large_factor: SCS failed to produce SCS_SOLVED",
            exitflag == SCS_SOLVED);

  /* ||L x - b||_inf, relative to ||b||_inf */
  res = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    res[j] = -d->b[j];
  }
  SCS(accum_by_a)(A, sol->x, res);
  for (j = 0; j < n; ++j) {
    res_max = MAX(res_max, ABS(res[j]));
  }
  scs_printf("small_kkt_large_factor: ||Lx - b||_inf = %.3e\n", res_max);
  scs_free(res);
  if (res_max > 1e-4 * b_max) {
    fail = "small_kkt_large_factor: ||L x - b|| too large";
  }
  if (!fail) {
    fail = verify_solution_correct(d, cone, stgs, &info, sol, exitflag);
  }

  SCS(free_data)(d);
  SCS(free_cone)(cone);
  SCS(free_sol)(sol);
  scs_free(stgs);
  return fail;
#endif
}
