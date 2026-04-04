#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * SOC size-1 test: proj_soc with q=1 reduces to x[0] = MAX(x[0], 0), i.e.
 * the non-negative orthant for a single scalar.
 * Problem: min t  s.t.  t in SOC_1 (t >= 0),  t = 1.5
 * Optimal: t = 1.5, obj = 1.5.
 *
 * Variables [t] (n=1), m=2:
 *   Row 0 (zero): t = 1.5  -> A[0,0]=1, b[0]=1.5
 *   Row 1 (SOC_1): s[1] = t  -> A[1,0]=-1, b[1]=0
 */
static const char *test_soc_size1(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 1.5;

  scs_float Ax[] = {1.0, -1.0};
  scs_int Ai[]   = {0, 1};
  scs_int Ap[]   = {0, 2};

  scs_float b[] = {1.5, 0.0};
  scs_float c[] = {1.0};

  scs_int q[] = {1};
  scs_int m = 2, n = 1;

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
  k->qsize = 1;
  k->q = q;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_soc_size1: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

/*
 * SOC size-2 test: exercises the q=2 fast path in proj_soc.
 * min t  s.t. (t, x) in SOC_2 (t >= |x|),  x = 1.5
 * Optimal: t = 1.5, obj = 1.5.
 *
 * Variables [t, x] (n=2), m=3:
 *   Row 0 (zero): x = 1.5  → A[0,1]=1, b[0]=1.5
 *   Rows 1-2 (SOC_2): s[1:3] = (t, x) in SOC_2
 *     row 1: -t+s1=0  → s1=t
 *     row 2: -x+s2=0  → s2=x
 */
static const char *test_soc_size2(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 1.5;

  scs_float Ax[] = {-1.0, 1.0, -1.0};
  scs_int Ai[]   = {1, 0, 2};
  scs_int Ap[]   = {0, 1, 3};

  scs_float b[] = {1.5, 0.0, 0.0};
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

  k->z = 1;
  k->qsize = 1;
  k->q = q;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_soc_size2: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

/*
 * SOC size-3 test: exercises the q=3 fast path in proj_soc.
 * min t  s.t. (t, x1, x2) in SOC_3  (t >= sqrt(x1^2+x2^2)),
 *              x1 = 1,  x2 = 2
 * Optimal: t = sqrt(5), obj = sqrt(5).
 *
 * Variables [t, x1, x2] (n=3), m=5:
 *   Row 0 (zero): x1 = 1   → A[0,1]=1, b[0]=1
 *   Row 1 (zero): x2 = 2   → A[1,2]=1, b[1]=2
 *   Rows 2-4 (SOC_3): s[2:5] = (t,x1,x2) in SOC_3
 *     row 2: -t+s2=0, row 3: -x1+s3=0, row 4: -x2+s4=0
 */
static const char *test_soc_size3(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 2.2360679774997896; /* sqrt(5) */

  scs_float Ax[] = {-1.0, 1.0, -1.0, 1.0, -1.0};
  scs_int Ai[]   = {2, 0, 3, 1, 4};
  scs_int Ap[]   = {0, 1, 3, 5};

  scs_float b[] = {1.0, 2.0, 0.0, 0.0, 0.0};
  scs_float c[] = {1.0, 0.0, 0.0};

  scs_int q[] = {3};
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
  k->qsize = 1;
  k->q = q;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_soc_size3: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}

/*
 * SOC size-5 test: exercises the general (q>3) code path in proj_soc.
 * min t  s.t.  (t, x1, x2, x3, x4) in SOC_5  (t >= ||(x1..x4)||),
 *              xi = 1 for i=1..4.
 * Optimal: t = sqrt(4) = 2, obj = 2.
 *
 * Variables [t, x1, x2, x3, x4] (n=5), m=9:
 *   Rows 0-3 (zero, 4 equalities): xi = 1
 *   Rows 4-8 (SOC_5): (t, x1, ..., x4) in SOC_5
 *
 * A cols (CSC):
 *   col 0 (t):  row 4 -> -1
 *   col 1 (x1): row 0 -> +1, row 5 -> -1
 *   col 2 (x2): row 1 -> +1, row 6 -> -1
 *   col 3 (x3): row 2 -> +1, row 7 -> -1
 *   col 4 (x4): row 3 -> +1, row 8 -> -1
 */
static const char *test_soc_size5(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 2.0; /* sqrt(1^2 + 1^2 + 1^2 + 1^2) = 2 */

  scs_float Ax5[] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  scs_int Ai5[]   = {4, 0, 5, 1, 6, 2, 7, 3, 8};
  scs_int Ap5[]   = {0, 1, 3, 5, 7, 9};

  scs_float b5[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  scs_float c5[] = {1.0, 0.0, 0.0, 0.0, 0.0};

  scs_int q5[] = {5};
  scs_int m = 9, n = 5;

  d->m = m;
  d->n = n;
  d->b = b5;
  d->c = c5;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = n;
  d->A->x = Ax5;
  d->A->i = Ai5;
  d->A->p = Ap5;

  k->z = 4;
  k->qsize = 1;
  k->q = q5;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_soc_size5: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
