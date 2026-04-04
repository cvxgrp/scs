#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Unbounded SOCP: min -(t + x)  s.t.  (t, x) in SOC_2.
 * On the boundary t = x >= 0, objective = -2t -> -inf.
 *
 * Variables [t, x] (n=2), m=2 (SOC_2):
 *   row 0: -t + s0 = 0  => s0 = t
 *   row 1: -x + s1 = 0  => s1 = x
 *
 * A cols (CSC):
 *   col 0 (t): row 0 -> -1
 *   col 1 (x): row 1 -> -1
 */
static const char *unbounded_socp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  scs_float Ax[] = {-1.0, -1.0};
  scs_int Ai[]   = {0, 1};
  scs_int Ap[]   = {0, 1, 2};

  scs_float b[] = {0.0, 0.0};
  scs_float c[] = {-1.0, -1.0};

  scs_int q[] = {2};
  scs_int m = 2, n = 2;

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

  k->qsize = 1;
  k->q = q;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("unbounded_socp: expected SCS_UNBOUNDED",
            exitflag == SCS_UNBOUNDED);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
