#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Power cone test: maximize t = (x1 * x2)^0.5 subject to x1 + x2 = 4.
 * Equivalently: min -t s.t. (x1, x2, t) in K_pow(0.5), x1+x2=4.
 * Optimal: x1 = x2 = 2, t = 2, obj = -2.
 *
 * SCS form: min c'v, Av + s = b, s in K
 * Variables v = [x1, x2, t] (n=3), m=4
 * Constraints:
 *   Row 0 (zero cone):  x1 + x2 = 4
 *   Rows 1-3 (power):   s[1:4] = (x1, x2, t) in K_pow(0.5)
 *
 * A cols (CSC):
 *   col 0 (x1): row 0 -> +1, row 1 -> -1
 *   col 1 (x2): row 0 -> +1, row 2 -> -1
 *   col 2 (t):  row 3 -> -1
 */
static const char *test_power_cone(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = -2.0;

  scs_float Ax[] = {1.0, -1.0, 1.0, -1.0, -1.0};
  scs_int Ai[]   = {0, 1, 0, 2, 3};
  scs_int Ap[]   = {0, 2, 4, 5};

  scs_float b[] = {4.0, 0.0, 0.0, 0.0};
  scs_float c[] = {0.0, 0.0, -1.0};

  scs_float p[] = {0.5};

  scs_int m = 4, n = 3;

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

  k->z = 1;
  k->psize = 1;
  k->p = p;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_power_cone: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  /* warm start */
  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("test_power_cone: warm-start failure", info.iter <= 25);
  success = ABS(info.pobj - opt) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_power_cone: SCS warm-start failed", success);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

/*
 * Infeasible power cone: (x, y, z) in K_pow(0.5) with x = -1 (impossible
 * since x must be >= 0 in the power cone).
 *
 * Variables v = [x, y, z] (n=3), m=4
 *   Row 0 (zero cone): x = -1   → A[0,0]=1, b[0]=-1, s[0]=0 → x=-1
 *   Rows 1-3 (power):  s[1:4] = (x, y, z) in K_pow(0.5)
 *
 * A:
 *   col 0 (x): row 0 -> +1, row 1 -> -1
 *   col 1 (y): row 2 -> -1
 *   col 2 (z): row 3 -> -1
 */
static const char *test_power_cone_infeasible(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  scs_float Ax[] = {1.0, -1.0, -1.0, -1.0};
  scs_int Ai[]   = {0, 1, 2, 3};
  scs_int Ap[]   = {0, 2, 3, 4};

  scs_float b[] = {-1.0, 0.0, 0.0, 0.0};
  scs_float c[] = {0.0, 0.0, -1.0};

  scs_float p[] = {0.5};

  scs_int m = 4, n = 3;

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

  k->z = 1;
  k->psize = 1;
  k->p = p;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_power_cone_infeasible: expected SCS_INFEASIBLE",
            exitflag == SCS_INFEASIBLE);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
