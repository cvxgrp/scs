#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Test box cone projection on a simple LP (no P matrix).
 *
 * K_box = {(t, s) : bl_i * t <= s_i <= bu_i * t, t >= 0}
 *
 * Problem: minimize x1 + x2
 *          subject to -1 <= x1 <= 1
 *                     -2 <= x2 <= 2
 * Solution: x1 = -1, x2 = -2, obj = -3.
 *
 * SCS formulation: Ax + s = b, s in K_box
 *   k->bsize = 3  (t + x1 + x2)
 *   k->bl = [-1, -2], k->bu = [1, 2]
 *   Row 0 (t-slot):  0*x1 + 0*x2 + s[0] = 1  =>  t = 1
 *   Row 1 (x1-slot): -x1 + s[1] = 0           =>  s[1] = x1
 *   Row 2 (x2-slot): -x2 + s[2] = 0           =>  s[2] = x2
 *
 * A (CSC, 3x2):
 *   col 0 (x1): row 1 -> -1
 *   col 1 (x2): row 2 -> -1
 */
static const char *test_box_cone_lp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = -3.0;

  scs_float Ax[] = {-1.0, -1.0};
  scs_int Ai[]   = {1, 2};
  scs_int Ap[]   = {0, 1, 2};

  scs_float b[] = {1.0, 0.0, 0.0};
  scs_float c[] = {1.0, 1.0};

  scs_float bl[] = {-1.0, -2.0};
  scs_float bu[] = {1.0, 2.0};

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

  k->bsize = 3;
  k->bl = bl;
  k->bu = bu;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_box_cone_lp: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
