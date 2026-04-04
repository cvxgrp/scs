#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Unbounded LP: min -x  s.t. x >= 0.
 * x is unbounded above, so the objective -> -inf.
 *
 * SCS form: Ax + s = b, s in K_nn (l=1)
 *   row 0: -x + s = 0, s >= 0  =>  s = x >= 0
 *
 * A = [[-1]], b = [0], c = [-1], k->l = 1.
 */
static const char *unbounded_lp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  scs_float Ax[] = {-1.0};
  scs_int Ai[]   = {0};
  scs_int Ap[]   = {0, 1};

  scs_float b[] = {0.0};
  scs_float c[] = {-1.0};

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

  k->l = m;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("unbounded_lp: expected SCS_UNBOUNDED",
            exitflag == SCS_UNBOUNDED);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
