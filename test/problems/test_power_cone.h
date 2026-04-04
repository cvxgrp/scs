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
 * Power cone test with p=0.9: exercises a non-symmetric power parameter.
 * Problem: maximize t = x^0.9 * y^0.1  s.t.  x = 1,  y = 1.
 * K_pow(0.9): x^0.9 * y^0.1 >= |t|, x,y >= 0.
 * Optimal: t* = 1^0.9 * 1^0.1 = 1,  obj = -1.
 *
 * Variables [x, y, t] (n=3), m=5:
 *   Row 0 (zero): x = 1  -> A[0,0]=1, b[0]=1
 *   Row 1 (zero): y = 1  -> A[1,1]=1, b[1]=1
 *   Rows 2-4 (power p=0.9): (x, y, t) in K_pow(0.9)
 *
 * A cols (CSC):
 *   col 0 (x): row 0 -> +1, row 2 -> -1
 *   col 1 (y): row 1 -> +1, row 3 -> -1
 *   col 2 (t): row 4 -> -1
 */
static const char *test_power_cone_p09(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = -1.0; /* 1^0.9 * 1^0.1 = 1 */

  scs_float Ax_p[] = {1.0, -1.0, 1.0, -1.0, -1.0};
  scs_int Ai_p[]   = {0, 2, 1, 3, 4};
  scs_int Ap_p[]   = {0, 2, 4, 5};

  scs_float b_p[] = {1.0, 1.0, 0.0, 0.0, 0.0};
  scs_float c_p[] = {0.0, 0.0, -1.0};

  scs_float pp[] = {0.9};

  scs_int m = 5, n = 3;

  d->m = m;
  d->n = n;
  d->b = b_p;
  d->c = c_p;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = n;
  d->A->x = Ax_p;
  d->A->i = Ai_p;
  d->A->p = Ap_p;

  k->z = 2;
  k->psize = 1;
  k->p = pp;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_power_cone_p09: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

/*
 * Dual power cone test (p < 0 in SCS signals dual cone).
 *
 * For p[i] < 0, cones.c projects onto the dual of K_pow(-p[i]) via Moreau:
 * Pi_{K*}(x) = x + Pi_K(-x).
 *
 * K_pow(0.5)* = {(u,v,w): 2*sqrt(u*v) >= |w|, u,v >= 0}.
 * Problem: maximize w  s.t.  (u, v, w) in K_pow(0.5)*,  u = 1,  v = 1.
 * With u=v=1: max w = 2*sqrt(1) = 2,  obj = -2.
 *
 * Variables [u, v, w] (n=3), m=5:
 *   Row 0 (zero): u = 1  -> A[0,0]=1, b[0]=1
 *   Row 1 (zero): v = 1  -> A[1,1]=1, b[1]=1
 *   Rows 2-4 (dual power p=-0.5): (u, v, w) in K_pow(0.5)*
 *
 * A cols (CSC):
 *   col 0 (u): row 0 -> +1, row 2 -> -1
 *   col 1 (v): row 1 -> +1, row 3 -> -1
 *   col 2 (w): row 4 -> -1
 */
static const char *test_dual_power_cone(void) {
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

  scs_float Ax_d[] = {1.0, -1.0, 1.0, -1.0, -1.0};
  scs_int Ai_d[]   = {0, 2, 1, 3, 4};
  scs_int Ap_d[]   = {0, 2, 4, 5};

  scs_float b_d[] = {1.0, 1.0, 0.0, 0.0, 0.0};
  scs_float c_d[] = {0.0, 0.0, -1.0};

  scs_float pd[] = {-0.5}; /* negative = dual of K_pow(0.5) */

  scs_int m = 5, n = 3;

  d->m = m;
  d->n = n;
  d->b = b_d;
  d->c = c_d;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = n;
  d->A->x = Ax_d;
  d->A->i = Ai_d;
  d->A->p = Ap_d;

  k->z = 2;
  k->psize = 1;
  k->p = pd;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_dual_power_cone: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

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

/*
 * Multiple power cones test (psize=2): exercises the loop over k->psize in
 * proj_cones, including both primal and dual cone variants.
 *
 * Two independent power cone constraints:
 *   Cone 1 (p=0.5):  (x1, x2, t1) in K_pow(0.5),  x1=1, x2=1 -> t1*=1
 *   Cone 2 (p=-0.5): (y1, y2, t2) in K_pow(0.5)*, y1=1, y2=1 -> t2*=2
 * (Dual K_pow(0.5)*: 2*sqrt(y1*y2) >= |t2|, so t2* = 2*sqrt(1) = 2.)
 *
 * Problem: minimize -(t1 + t2)
 * Optimal: t1* = 1, t2* = 2, obj* = -3.
 *
 * Variables [x1, x2, t1, y1, y2, t2] (n=6), m=10:
 *   Rows 0-3 (z=4):  x1=1, x2=1, y1=1, y2=1
 *   Rows 4-6 (p=0.5):  (x1, x2, t1)
 *   Rows 7-9 (p=-0.5): (y1, y2, t2)
 *
 * A cols (CSC):
 *   col 0 (x1): row 0 -> +1, row 4 -> -1
 *   col 1 (x2): row 1 -> +1, row 5 -> -1
 *   col 2 (t1): row 6 -> -1
 *   col 3 (y1): row 2 -> +1, row 7 -> -1
 *   col 4 (y2): row 3 -> +1, row 8 -> -1
 *   col 5 (t2): row 9 -> -1
 */
static const char *test_multi_power(void) {
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

  scs_float Axmp[] = {1.0,-1.0, 1.0,-1.0, -1.0, 1.0,-1.0, 1.0,-1.0, -1.0};
  scs_int Aimp[]   = {0, 4, 1, 5, 6, 2, 7, 3, 8, 9};
  scs_int Apmp[]   = {0, 2, 4, 5, 7, 9, 10};

  scs_float bmp[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  scs_float cmp[] = {0.0, 0.0, -1.0, 0.0, 0.0, -1.0};

  scs_float pmp[] = {0.5, -0.5}; /* primal K_pow(0.5) then dual K_pow(0.5)* */

  scs_int m = 10, n = 6;

  d->m = m;
  d->n = n;
  d->b = bmp;
  d->c = cmp;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = n;
  d->A->x = Axmp;
  d->A->i = Aimp;
  d->A->p = Apmp;

  k->z = 4;
  k->psize = 2;
  k->p = pmp;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_multi_power: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
