#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Test scs_update (b and/or c) for a simple LP.
 *
 * Base problem: min x  s.t.  x >= 2
 *   row 0: -x + s = -2, s >= 0  =>  x >= 2
 *   n=1, m=1, k->l=1
 *   A=[[-1]], b=[-2], c=[1]
 *   Optimal: x=2, obj=2.
 *
 * After scs_update(b=[-3]): x >= 3, obj=3.
 * After scs_update(c=[2]):  min 2x s.t. x >= 3, obj=6.
 * After scs_update(b=[-2], c=[1]): back to original, obj=2.
 */
static const char *lp_update(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_int success;
  const char *fail;

  scs_float Ax[] = {-1.0};
  scs_int Ai[]   = {0};
  scs_int Ap[]   = {0, 1};

  scs_float b[] = {-2.0};
  scs_float c[] = {1.0};

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

  ScsWork *w = scs_init(d, k, stgs);
  mu_assert("lp_update: scs_init returned NULL", w != SCS_NULL);

  /* base solve: min x s.t. x >= 2, expected obj=2 */
  exitflag = scs_solve(w, sol, &info, 0);
  success = exitflag == SCS_SOLVED && ABS(info.pobj - 2.0) < 1e-4;
  mu_assert("lp_update: base solve failed", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) return fail;

  /* update b: x >= 3, expected obj=3 */
  b[0] = -3.0;
  scs_update(w, b, SCS_NULL);
  exitflag = scs_solve(w, sol, &info, 1);
  success = exitflag == SCS_SOLVED && ABS(info.pobj - 3.0) < 1e-4;
  mu_assert("lp_update: b-update solve failed", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) return fail;

  /* update c: min 2x s.t. x >= 3, expected obj=6 */
  c[0] = 2.0;
  scs_update(w, SCS_NULL, c);
  exitflag = scs_solve(w, sol, &info, 1);
  success = exitflag == SCS_SOLVED && ABS(info.pobj - 6.0) < 1e-4;
  mu_assert("lp_update: c-update solve failed", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) return fail;

  /* update both b and c back to original: min x s.t. x >= 2, expected obj=2 */
  b[0] = -2.0;
  c[0] = 1.0;
  scs_update(w, b, c);
  exitflag = scs_solve(w, sol, &info, 1);
  success = exitflag == SCS_SOLVED && ABS(info.pobj - 2.0) < 1e-4;
  mu_assert("lp_update: b+c-update solve failed", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) return fail;

  /* null/null update: no change to problem, result must be identical */
  scs_update(w, SCS_NULL, SCS_NULL);
  exitflag = scs_solve(w, sol, &info, 1);
  success = exitflag == SCS_SOLVED && ABS(info.pobj - 2.0) < 1e-4;
  mu_assert("lp_update: null-null update changed result", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  scs_finish(w);
  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
