#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * The box cone's first entry t is weighted like an equality row: its R_y entry
 * is 1 / (BOX_T_RY_FACTOR * scale), the zero cone's weight, where every other
 * row outside the zero cone gets 1 / scale (SCS(set_r_y), glbopts.h). With
 * adaptive_diag_scale each row's entry is further divided by its own
 * multiplier.
 *
 * In the box cone's usual form, l <= a'x <= u written as (1, a'x) in
 * {(t, s) : t l <= s <= t u}, the row of A for t is zero and its entry of b
 * is 1, so the row says t = 1. Weighted like an inequality row, t could stop
 * off 1 by about SCS's primal tolerance, eps (1 + ||A x||), and every bound is
 * then met only up to |t - 1| times the bound.
 *
 * The problems are built from HS76 (Hock and Schittkowski, test problem 76),
 * a convex QP in 4 variables:
 *
 *   minimize    x1^2 + x2^2 / 2 + x3^2 + x4^2 / 2 - x1 x3 + x3 x4
 *               - x1 - 3 x2 + x3 - x4
 *   subject to  x1 + 2 x2 + x3 + x4 <= 5
 *               3 x1 + x2 + 2 x3 - x4 <= 4
 *               x2 + 4 x3 >= 1.5
 *               x >= 0,
 *
 * with optimal value -103 / 22 = -4.6818... at x = (3, 23, 0, 6) / 11, where
 * the first constraint and x3 >= 0 are the active ones.
 */

#define HS76_OPT ((scs_float)(-103. / 22.))

/* the three forms of HS76 below, m = 8 rows each */
enum hs76_form {
  /* as a QP benchmark harness poses it: every row in one box cone,
   * bsize = 8. Row 0 is t (A row zero, b = 1), rows 1 to 7 hold a'x for the
   * seven constraints, with their bounds in bl, bu. */
  HS76_BOX,
  /* the first constraint, active at the optimum, as an equality (row 0, a
   * zero-cone row), then a box cone of size 7: t is row 1, rows 2 to 7 hold
   * the other six constraints. */
  HS76_ZERO_BOX,
  /* no box cone with bounds: the first constraint as an equality, the other
   * six as nonnegative rows, and row 7 a box cone of size 1 (t >= 0 alone,
   * A row zero, b = 1), which has no bounds and is weighted like a
   * nonnegative row. */
  HS76_NO_BOX
};

/* P (upper triangle) and c of HS76 */
static scs_float hs76_Px[] = {2., 1., -1., 2., 1., 1.};
static scs_int hs76_Pi[] = {0, 1, 0, 2, 2, 3};
static scs_int hs76_Pp[] = {0, 1, 2, 4, 6};
static scs_float hs76_c[] = {-1., -3., 1., -1.};
/* A has 14 nonzeros in every form, 3, 4, 4 and 3 per column */
static scs_int hs76_Ap[] = {0, 3, 7, 11, 14};

/* A x + s = b, s in K. In the box forms a box row reads s_j = a_j'x (row of A
 * -a_j, b_j = 0) and the equality row a_1'x = 5 (row of A a_1, b = 5). */
static scs_float hs76_box_Ax[] = {-1., -3., -1., -2., -1., -1., -1.,
                                  -1., -2., -4., -1., -1., 1.,  -1.};
static scs_int hs76_box_Ai[] = {1, 2, 4, 1, 2, 3, 5, 1, 2, 3, 6, 1, 2, 7};
static scs_float hs76_box_b[] = {1., 0., 0., 0., 0., 0., 0., 0.};
static scs_float hs76_box_bl[] = {-INFINITY, -INFINITY, 1.5, 0., 0., 0., 0.};
static scs_float hs76_box_bu[] = {5.,       4.,       INFINITY, INFINITY,
                                  INFINITY, INFINITY, INFINITY};

static scs_float hs76_zero_box_Ax[] = {1., -3., -1., 2., -1., -1., -1.,
                                       1., -2., -4., -1., 1., 1.,  -1.};
static scs_int hs76_zero_box_Ai[] = {0, 2, 4, 0, 2, 3, 5, 0, 2, 3, 6, 0, 2, 7};
static scs_float hs76_zero_box_b[] = {5., 1., 0., 0., 0., 0., 0., 0.};
static scs_float hs76_zero_box_bl[] = {-INFINITY, 1.5, 0., 0., 0., 0.};
static scs_float hs76_zero_box_bu[] = {4.,       INFINITY, INFINITY,
                                       INFINITY, INFINITY, INFINITY};

/* nonnegative rows s = b - A x >= 0: 4 - (3 x1 + x2 + 2 x3 - x4),
 * x2 + 4 x3 - 1.5 and x itself; row 7 is s = 1 */
static scs_float hs76_no_box_Ax[] = {1., 3., -1., 2., 1., -1., -1.,
                                     1., 2., -4., -1., 1., -1., -1.};
static scs_int hs76_no_box_Ai[] = {0, 1, 3, 0, 1, 2, 4, 0, 1, 2, 5, 0, 1, 6};
static scs_float hs76_no_box_b[] = {5., 4., -1.5, 0., 0., 0., 0., 1.};

/* Point d and k at HS76 in the given form; d->A comes from scs_calloc and the
 * caller frees it, everything else is static. */
static void hs76_set(ScsData *d, ScsCone *k, enum hs76_form form) {
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsMatrix *P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  memset(k, 0, sizeof(*k));
  A->m = 8;
  A->n = 4;
  A->p = hs76_Ap;
  P->m = 4;
  P->n = 4;
  P->x = hs76_Px;
  P->i = hs76_Pi;
  P->p = hs76_Pp;
  d->m = 8;
  d->n = 4;
  d->A = A;
  d->P = P;
  d->c = hs76_c;
  switch (form) {
  case HS76_BOX:
    A->x = hs76_box_Ax;
    A->i = hs76_box_Ai;
    d->b = hs76_box_b;
    k->bsize = 8;
    k->bl = hs76_box_bl;
    k->bu = hs76_box_bu;
    break;
  case HS76_ZERO_BOX:
    A->x = hs76_zero_box_Ax;
    A->i = hs76_zero_box_Ai;
    d->b = hs76_zero_box_b;
    k->z = 1;
    k->bsize = 7;
    k->bl = hs76_zero_box_bl;
    k->bu = hs76_zero_box_bu;
    break;
  case HS76_NO_BOX:
    A->x = hs76_no_box_Ax;
    A->i = hs76_no_box_Ai;
    d->b = hs76_no_box_b;
    k->z = 1;
    k->l = 6;
    k->bsize = 1;
    break;
  }
}

/* Check every R_y entry of the workspace against its weight: rho_y_i times
 * scale times row i's multiplier must be 1 / ZERO_CONE_RY_FACTOR for a
 * zero-cone row, 1 / BOX_T_RY_FACTOR for the t entry of a box cone of size 2
 * or more, and 1 for every other row. A wrong weight is off by a factor of
 * about 1000, so the relative tolerance is loose. */
static const char *check_r_y_weights(const ScsWork *w, const ScsCone *k,
                                     const char *when) {
  scs_int i, n = w->d->n, m = w->d->m;
  scs_int t_row = (k->bsize > 1) ? k->z + k->l : -1;
  for (i = 0; i < m; ++i) {
    scs_float mult = w->scale_mults ? w->scale_mults[i] : 1.;
    scs_float base = w->diag_r[n + i] * w->stgs->scale * mult;
    scs_float want = 1., one = 1.;
    if (i < k->z) {
      want = (scs_float)(1. / ZERO_CONE_RY_FACTOR);
    } else if (i == t_row) {
      want = (scs_float)(1. / BOX_T_RY_FACTOR);
    }
    if (ABS(base / want - one) > 1e-3) {
      scs_printf("%s: row %li: rho_y * scale * multiplier = %.6e, expected "
                 "%.6e (scale %.3e, multiplier %.3e)\n",
                 when, (long)i, base, want, w->stgs->scale, mult);
      return "R_y entry does not match its cone's weight";
    }
  }
  return SCS_NULL;
}

/* R_y of each form at setup and after scale updates. The solve starts from
 * scale = 1e-5, three orders of magnitude below where this problem settles
 * (about 1e-2), so the first update (after RESCALING_MIN_ITERS iterations) is
 * certain; eps is set beyond reach after setup (in single precision too) so
 * that the solve runs long enough to make it. */
static const char *r_y_weights_one(enum hs76_form form, const char *name) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  ScsWork *w;
  const char *fail = SCS_NULL;

  hs76_set(d, k, form);
  scs_set_default_settings(stgs);
  stgs->verbose = 0;
  stgs->scale = 1e-5;
  stgs->adaptive_diag_scale = 1;

  w = scs_init(d, k, stgs);
  if (!w) {
    fail = "scs_init returned NULL";
  } else {
    fail = check_r_y_weights(w, k, name);
    if (!fail) {
      w->stgs->eps_abs = 1e-12;
      w->stgs->eps_rel = 1e-12;
      w->stgs->max_iters = 400;
      scs_solve(w, sol, &info, 0);
      scs_printf("%s: %li iterations, %li scale updates, final scale %.3e\n",
                 name, (long)info.iter, (long)info.scale_updates,
                 w->stgs->scale);
      if (info.scale_updates < 1) {
        fail = "no scale update, the check after updates is vacuous";
      } else {
        fail = check_r_y_weights(w, k, name);
      }
    }
    scs_finish(w);
  }

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(d);
  scs_free(k);
  scs_free(stgs);
  return fail;
}

/* The box cone's t entry gets the zero cone's weight, at setup and after
 * scale updates, and every other row keeps its own. */
static const char *test_box_t_weight_r_y(void) {
  return r_y_weights_one(HS76_ZERO_BOX, "test_box_t_weight_r_y");
}

/* Without a box cone with bounds (a box cone of size 1 included) every R_y
 * entry is what it was before the box t entry got its weight. */
static const char *test_box_t_weight_no_box(void) {
  return r_y_weights_one(HS76_NO_BOX, "test_box_t_weight_no_box");
}

/* HS76 as a benchmark suite poses it, every row in a box cone, at the default
 * eps of 1e-4. Weighted like an inequality row, t stopped at 1 + 4.1e-4 here
 * (macOS, Accelerate and QDLDL, any starting scale near the default), close to
 * SCS's primal tolerance of about 1e-4 (1 + 5), and the first constraint was
 * violated by 2.2e-3; weighted like an equality row, t stops within 1e-9 of 1
 * and the violation is below 1e-6. The test asks for |t - 1| below 1e-5,
 * between the two by orders of magnitude, so it fails if t is weighted like an
 * inequality row again but not on the roundoff of another platform or of
 * single precision. */
static const char *test_box_t_weight_stop(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag, j;
  scs_float t_err, viol = 0., one = 1.;
  const char *fail;
  /* the seven constraints a'x in [l, u] of HS76 */
  static const scs_float a[7][4] = {{1., 2., 1., 1.}, {3., 1., 2., -1.},
                                    {0., 1., 4., 0.}, {1., 0., 0., 0.},
                                    {0., 1., 0., 0.}, {0., 0., 1., 0.},
                                    {0., 0., 0., 1.}};

  hs76_set(d, k, HS76_BOX);
  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-4;
  stgs->eps_rel = 1e-4;
  stgs->verbose = 0;

  exitflag = scs(d, k, stgs, sol, &info);

  /* t is s[0]; s is returned with tau divided out */
  t_err = ABS(sol->s[0] - one);
  for (j = 0; j < 7; ++j) {
    scs_float ax = a[j][0] * sol->x[0] + a[j][1] * sol->x[1] +
                   a[j][2] * sol->x[2] + a[j][3] * sol->x[3];
    viol = MAX(viol, MAX(hs76_box_bl[j] - ax, 0.) + MAX(ax - hs76_box_bu[j], 0.));
  }
  scs_printf("test_box_t_weight_stop: %li iterations, |t - 1| = %.3e, largest "
             "bound violation %.3e, objective error %.3e\n",
             (long)info.iter, t_err, viol, ABS(info.pobj - HS76_OPT));

  mu_assert("test_box_t_weight_stop: expected SCS_SOLVED",
            exitflag == SCS_SOLVED);
  mu_assert_less("test_box_t_weight_stop: |t - 1| at the stop", t_err, 1e-5);
  mu_assert_less("test_box_t_weight_stop: objective error",
                 ABS(info.pobj - HS76_OPT), SOL_TOL);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(d);
  scs_free(k);
  scs_free(stgs);
  return fail;
}
