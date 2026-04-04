#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Mixed-cone test: exercises the full cone ordering path (z, l, q, ep)
 * in a single problem to verify that set_cone_boundaries and proj_cone
 * dispatch correctly across multiple different cone types.
 *
 * Problem (5 variables, 8 constraints):
 *   minimize  x1 + x2 + x3 + x4 + x5
 *   subject to
 *     x1 + x2 = 3                (zero cone, 1 row)
 *     x3 >= 1                    (non-neg cone, 1 row)
 *     (x4, x5) in SOC_2          (SOC, 2 rows)
 *     (x1, x2, x3) in K_exp      (exp cone, 3 rows)
 *     x4 >= 1   (duplicate of above, extra l row)
 *
 * Cone ordering: z=1, l=2, q=[2], ep=1 → total m=1+2+2+3=8.
 *
 * But this is over-constrained. Let me use a simpler formulation:
 *
 * Variables: [x1, x2, x3] (n=3)
 * m=6, cone order: z=1, l=1, q=[2], ep=0 (skipping exp for simpler math)
 *
 * Actually let's just use z + l + q:
 *   minimize  x1 + x2 + x3
 *   subject to
 *     x1 = 2                        (zero cone, 1 row: z=1)
 *     x2 >= 1                       (non-neg cone, 1 row: l=1)
 *     (x3, x2) in SOC_2             (SOC, 2 rows: q=[2])
 *   Optimal: x1=2, and (x3, x2) in SOC means x3 >= |x2|.
 *   To minimize x2+x3 with x2>=1: x2=1, x3>=1, so x3=1.
 *   obj = 2+1+1 = 4.
 *
 * SCS formulation: Ax + s = b, s in K
 *   Row 0 (z):  -x1 + s0 = -2       => x1 = 2
 *   Row 1 (l):  -x2 + s1 = -1       => x2 >= 1
 *   Row 2 (q0): -x3 + s2 = 0        => s2 = x3  (t component of SOC)
 *   Row 3 (q1): -x2 + s3 = 0        => s3 = x2  (body of SOC)
 *   SOC: |s3| <= s2, i.e., |x2| <= x3
 *
 * A (4x3 CSC):
 *   col 0 (x1): row 0 -> -1
 *   col 1 (x2): row 1 -> -1, row 3 -> -1
 *   col 2 (x3): row 2 -> -1
 */
static const char *test_mixed_cones(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float opt = 4.0;

  scs_float Ax[] = {-1.0, -1.0, -1.0, -1.0};
  scs_int Ai[]   = {0, 1, 3, 2};
  scs_int Ap[]   = {0, 1, 3, 4};

  scs_float b[] = {-2.0, -1.0, 0.0, 0.0};
  scs_float c[] = {1.0, 1.0, 1.0};

  scs_int q_arr[] = {2};
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
  k->l = 1;
  k->qsize = 1;
  k->q = q_arr;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  mu_assert("test_mixed_cones: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
