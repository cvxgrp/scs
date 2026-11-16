#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * cvxgrp/scs#531: with Apple Accelerate, dsyevr reports no convergence
 * (info = 1) on the 11 x 11 block below, which 3.3.2 turned into a failure
 * of the whole solve ("error in project_cones" after iteration 0). The block
 * is the first PSD iterate of the problem in issue_531 (CVXPY's
 * matrix_frac(x, P) s.t. x == 5, P == eye(11)): close to a*I + b*11' with
 * a < 0, so ten of its eigenvalues agree to the last few bits, which is the
 * case MRRR struggles with. The exact bits matter (a clean a*I + b*11' does
 * not trigger it), so they are stored as hex floats. The projection now
 * falls back to syev when syevr does not converge.
 *
 * Other LAPACKs (OpenBLAS, netlib) converge here, so on those the tests only
 * check the result; they exercise the fallback where dsyevr fails.
 */

/* svec of the block, every eigenvalue is negative (about -0.1776 and
 * -0.1702), so its projection onto the PSD cone is zero */
static const double issue_531_block[] = {
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cedap-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cedfp-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cedfp-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cedep-11,
    0x1.63423add1ceep-11,
    0x1.63423add1ceep-11,
    0x1.63423add1cedfp-11,
    0x1.63423add1cedap-11,
    0x1.63423add1ceep-11,
    0x1.63423add1cedap-11,
    0x1.63423add1cedfp-11,
    0x1.63423add1ceep-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee8p-11,
    0x1.63423add1cee8p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee8p-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1ceep-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee3p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1ceep-11,
    0x1.63423add1ceep-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cedfp-11,
    0x1.63423add1cee3p-11,
    0x1.63423add1cee2p-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cee8p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee6p-11,
    0x1.63423add1cee2p-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cedfp-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1ceep-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cedep-11,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee8p-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1cee2p-11,
    0x1.63423add1cee5p-11,
    -0x1.0042a2d83de55p-3,
    0x1.63423add1ceep-11,
    -0x1.0042a2d83de55p-3};

static const char *issue_531_proj(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsConeWork *cw;
  scs_int s_dims[1] = {11};
  scs_float x[66];
  scs_int i, ret;

  k->ssize = 1;
  k->s = s_dims;
  cw = SCS(init_cone)(k, 66);
  mu_assert("issue_531_proj: cone init failed", cw);
  for (i = 0; i < 66; ++i) {
    x[i] = (scs_float)issue_531_block[i];
  }
  ret = SCS(proj_dual_cone)(x, cw, SCS_NULL, SCS_NULL);
  SCS(finish_cone)(cw);
  scs_free(k);
  mu_assert("issue_531_proj: PSD projection failed", ret == 0);
  mu_assert("issue_531_proj: projection of a negative definite block is not 0",
            SCS(norm_inf)(x, 66) <= VERIFY_TOL(1e-12, 1e-5));
  return 0;
}

/* entries of A are 1, -1 and -sqrt(2) */
#define ISSUE_531_SQRT2 (1.4142135623730951)

/* The whole problem from the issue: m = 276, n = 78, 132 zero rows and PSD
 * cones of sizes 11 and 12, optimal value 275. */
static const char *issue_531(void) {
#ifdef SFLOAT
  /* The bug is in the projection, which issue_531_proj covers in single
   * precision too. The whole solve does not reach the default 1e-4 here in
   * float arithmetic: with the indirect solver it hits max_iters, and the
   * outcome depends on the platform's rounding. */
  return mu_skipped;
#endif
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag, i;
  const char *fail;

  scs_int Ap[] = {
      0, 1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 47, 51, 55, 59, 63,
      67, 71, 75, 79, 83, 86, 90, 94, 98, 102, 106, 110, 114, 118, 121, 125,
      129, 133, 137, 141, 145, 149, 152, 156, 160, 164, 168, 172, 176, 179,
      183, 187, 191, 195, 199, 202, 206, 210, 214, 218, 221, 225, 229, 233,
      236, 240, 244, 247, 251, 254, 256, 258, 260, 262, 264, 266, 268, 270,
      272, 274, 276};
  scs_int Ai[] = {
      275, 11, 132, 198, 12, 22, 133, 199, 13, 33, 134, 200, 14, 44, 135,
      201, 15, 55, 136, 202, 16, 66, 137, 203, 17, 77, 138, 204, 18, 88,
      139, 205, 19, 99, 140, 206, 20, 110, 141, 207, 21, 121, 142, 208, 23,
      143, 210, 24, 34, 144, 211, 25, 45, 145, 212, 26, 56, 146, 213, 27,
      67, 147, 214, 28, 78, 148, 215, 29, 89, 149, 216, 30, 100, 150, 217,
      31, 111, 151, 218, 32, 122, 152, 219, 35, 153, 221, 36, 46, 154, 222,
      37, 57, 155, 223, 38, 68, 156, 224, 39, 79, 157, 225, 40, 90, 158,
      226, 41, 101, 159, 227, 42, 112, 160, 228, 43, 123, 161, 229, 47, 162,
      231, 48, 58, 163, 232, 49, 69, 164, 233, 50, 80, 165, 234, 51, 91,
      166, 235, 52, 102, 167, 236, 53, 113, 168, 237, 54, 124, 169, 238, 59,
      170, 240, 60, 70, 171, 241, 61, 81, 172, 242, 62, 92, 173, 243, 63,
      103, 174, 244, 64, 114, 175, 245, 65, 125, 176, 246, 71, 177, 248, 72,
      82, 178, 249, 73, 93, 179, 250, 74, 104, 180, 251, 75, 115, 181, 252,
      76, 126, 182, 253, 83, 183, 255, 84, 94, 184, 256, 85, 105, 185, 257,
      86, 116, 186, 258, 87, 127, 187, 259, 95, 188, 261, 96, 106, 189, 262,
      97, 117, 190, 263, 98, 128, 191, 264, 107, 192, 266, 108, 118, 193,
      267, 109, 129, 194, 268, 119, 195, 270, 120, 130, 196, 271, 131, 197,
      273, 0, 209, 1, 220, 2, 230, 3, 239, 4, 247, 5, 254, 6, 260, 7, 265,
      8, 269, 9, 272, 10, 274};
  scs_float Ax[] = {
      -1., 1., -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2,
      -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1.,
      -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1.,
      -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1.,
      1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2,
      1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2,
      1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2,
      -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1.,
      -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1.,
      -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1.,
      -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1.,
      1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1.,
      1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2,
      1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2,
      1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2,
      1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2,
      1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1., 1., -ISSUE_531_SQRT2, -ISSUE_531_SQRT2, 1., -1., -1., 1.,
      -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2,
      1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2, 1., -ISSUE_531_SQRT2};
  scs_int b_idx[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 23, 35, 47, 59, 71, 83, 95, 107, 119, 131};
  scs_float b_val[] = {5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  scs_float b[276] = {0}, c[78] = {0};
  scs_int s[] = {11, 12};
  scs_int m = 276, n = 78;

  for (i = 0; i < (scs_int)(sizeof(b_idx) / sizeof(b_idx[0])); ++i) {
    b[b_idx[i]] = b_val[i];
  }
  c[0] = 1.;

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

  k->z = 132;
  k->s = s;
  k->ssize = 2;

  scs_set_default_settings(stgs);

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("issue_531: expected SCS_SOLVED", exitflag == SCS_SOLVED);
  mu_assert("issue_531: wrong optimal value",
            ABS(info.pobj - (scs_float)275.) <= (scs_float)2.75);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
