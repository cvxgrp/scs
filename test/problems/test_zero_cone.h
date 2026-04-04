#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Pure zero-cone problem: only equality constraints, no inequality cones.
 * The zero-cone projection (src/cones.c) simply memsets the slice to 0.
 *
 * Problem: minimize x1 + x2
 *          subject to x1 + x2 = 3
 *                     x1 - x2 = 1
 * Solution: x1 = 2, x2 = 1, obj = 3.
 *
 * Variables [x1, x2] (n=2), m=2, k->z=2 (only zero cone).
 *
 * A cols (CSC):
 *   col 0 (x1): row 0 -> +1, row 1 -> +1
 *   col 1 (x2): row 0 -> +1, row 1 -> -1
 */
static const char *test_zero_cone(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 3.0;

  scs_float Ax[] = {1.0, 1.0, 1.0, -1.0};
  scs_int Ai[]   = {0, 1, 0, 1};
  scs_int Ap[]   = {0, 2, 4};

  scs_float b[] = {3.0, 1.0};
  scs_float c[] = {1.0, 1.0};

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

  k->z = 2; /* only zero cone — no l, q, s, ep, ed, p */

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_zero_cone: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
