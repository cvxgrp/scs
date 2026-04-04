#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Dual exponential cone test.
 *
 * K_exp* (dual exp cone) = {(u,v,w): u<=0, -u*exp(v/u - 1) <= w}
 *                        U {u=0, v>=0, w>=0}
 * (Friberg 2021 convention used by SCS.)
 *
 * Problem: minimize w  s.t.  (u, v, w) in K_exp*,  u = -1,  v = 1.
 *
 * With u=-1, v=1 (first branch, u<0):
 *   -(-1)*exp(1/(-1) - 1) = exp(-2) <= w
 * Optimal: w* = exp(-2) ~= 0.13534, obj ~= 0.13534.
 *
 * Cone ordering: z=2 (equality rows 0-1), ed=1 (dual exp triple rows 2-4).
 * Variables [u, v, w] (n=3), m=5.
 *
 * Constraints:
 *   Row 0 (z): u = -1  -> A[0,0]=1, b[0]=-1
 *   Row 1 (z): v =  1  -> A[1,1]=1, b[1]= 1
 *   Row 2 (ed): s[2] = u  -> A[2,0]=-1, b[2]=0
 *   Row 3 (ed): s[3] = v  -> A[3,1]=-1, b[3]=0
 *   Row 4 (ed): s[4] = w  -> A[4,2]=-1, b[4]=0
 *
 * A cols (CSC):
 *   col 0 (u): row 0 -> +1, row 2 -> -1
 *   col 1 (v): row 1 -> +1, row 3 -> -1
 *   col 2 (w): row 4 -> -1
 */
static const char *test_dual_exp_cone(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  /* exp(-2) = 0.13533528323661270... */
  scs_float opt = 0.13533528323661270;

  scs_float Ax[] = {1.0, -1.0, 1.0, -1.0, -1.0};
  scs_int Ai[] = {0, 2, 1, 3, 4};
  scs_int Ap[] = {0, 2, 4, 5};

  scs_float b[] = {-1.0, 1.0, 0.0, 0.0, 0.0};
  scs_float c[] = {0.0, 0.0, 1.0};

  scs_int m = 5, n = 3;

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

  k->z = 2;
  k->ed = 1;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_dual_exp_cone: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
