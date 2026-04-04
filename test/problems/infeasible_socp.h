#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Infeasible SOCP: (t, x) in SOC_2  AND  t <= -1.
 * Since SOC_2 requires t >= |x| >= 0, t <= -1 is infeasible.
 *
 * Variables [t, x] (n=2), m=3.
 * Cone order: l=1 first, then q=[2].
 *   Row 0 (l=1):  t + s0 = -1, s0 >= 0  =>  t <= -1
 *   Row 1 (q[0]): -t + s1 = 0  => s1 = t
 *   Row 2 (q[1]): -x + s2 = 0  => s2 = x
 *   => (t, x) in SOC_2 AND t <= -1: infeasible since SOC_2 requires t >= 0.
 *
 * A cols (CSC):
 *   col 0 (t): row 0 -> +1, row 1 -> -1
 *   col 1 (x): row 2 -> -1
 */
static const char *infeasible_socp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  scs_float Ax[] = {1.0, -1.0, -1.0};
  scs_int Ai[]   = {0, 1, 2};
  scs_int Ap[]   = {0, 2, 3};

  scs_float b[] = {-1.0, 0.0, 0.0};
  scs_float c[] = {1.0, 0.0};

  scs_int q[] = {2};
  scs_int m = 3, n = 2;

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

  k->l = 1;
  k->qsize = 1;
  k->q = q;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("infeasible_socp: expected SCS_INFEASIBLE",
            exitflag == SCS_INFEASIBLE);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
