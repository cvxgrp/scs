#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Test PSD cone with n=1 (scalar SDP): exercises the fast path
 * proj_semi_definite_cone(X, n=1) which just does X[0] = MAX(X[0], 0)
 * without calling LAPACK eigendecomposition.  Works without USE_LAPACK.
 *
 * Problem: minimize x
 *          subject to x >= 1  (via 1x1 PSD cone)
 *
 * Formulation: Ax + s = b, s in S^1_+
 *   A = [[-1]], b = [-1], c = [1], n=1, m=1
 *   s = x - 1, s >= 0  =>  x >= 1
 *   Optimal: x = 1, obj = 1.
 *
 * k->ssize = 1, k->s = [1]   (one 1x1 PSD block)
 * PSD(1) contributes 1*(1+1)/2 = 1 entry in the constraint vector.
 */
static const char *test_psd_n1(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 1.0;

  scs_float Ax[] = {-1.0};
  scs_int Ai[]   = {0};
  scs_int Ap[]   = {0, 1};

  scs_float b[] = {-1.0};
  scs_float c[] = {1.0};

  scs_int s_arr[] = {1};   /* one 1x1 PSD block */
  scs_int m = 1, n = 1;

  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = n;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  k->ssize = 1;
  k->s = s_arr;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_psd_n1: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
