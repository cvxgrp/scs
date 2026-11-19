#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"
#include <string.h>

/*
 * The contract of scs_solve_lin_sys(w, b, s, tol, max_refine_steps, res)
 * (linsys.h) for the linear solver the test binary is built with, checked
 * against residuals the test computes itself from the problem data: the warm
 * start (direct solvers: the correction form and its fallback to a plain
 * solve; CG: its starting point), refinement (direct solvers: steps; CG:
 * restarts), and the residual res it reports, measured with the K of the
 * solver's rule. K = [R_x + P, A'; A, -R_y] with a random A (two thirds of
 * its entries, scaled so that ||K||_inf is about 10), a P with off-diagonal
 * entries and R in [0.5, 1.5], so that ||x|| is about ||b||, at two sizes:
 * about 8e4 nonzeros in K, from cuDSS's default CUDSS_MIN_PRODUCT_NNZ on (cuDSS
 * multiplies on the GPU there), and about 1e4, below it.
 *
 * Every tolerance is relative to ||b||_inf and lies between the rounding of
 * a residual (eps (||b|| + ||K|| ||x||): about 1e-14 of ||b|| in double
 * precision, 1e-6 in single) and what the bug it targets gives (1e-2 of
 * ||b|| or more: the residual of the wrong vector or with the wrong K, a
 * correction form that drops its warm start or uses another one), so both
 * precisions use the same constants. Each case prints its values.
 *
 * Phase 1, K as given to scs_init_lin_sys_work:
 *  (a) a plain solve: res is the test's ||b - K x|| and x is accurate;
 *  (b) from a warm start near the solution (1e-2 ||x|| away): res matches
 *      and x is accurate (a correction form that drops the warm start, or
 *      takes another vector for it, is off by about 1e-1 ||b||);
 *  (c) direct solvers, from a warm start 1e6 ||x|| away: the solver falls
 *      back to a plain solve (res at most 100 times the plain one's; the
 *      correction form from it would be about 1e6 times worse);
 *  (d) CG solvers, from a warm start so far (LSS_CG_FAR) that the rounding
 *      of b - K s leaves CG's tracked residual far below the true one (the
 *      premise, res > 100 tol, is asserted): two restarts bring res to at
 *      most 10 tol, and below the unrestarted res;
 *  (e) refinement (max_refine_steps 3; tol 0 for direct solvers, CG's tol
 *      for CG): never worse than the plain solve (2 res for run-to-run
 *      variation), and res matches.
 * Phase 2, after scs_update_lin_sys_diag_r with R' = R (1 + LSS_DELTA) in
 * another array: the solver factorizes (or, for CG, solves with) K' while
 * SCS's products keep measuring K, so a plain solve's residual with K is
 * about LSS_DELTA ||R x||. A solver that measures with its own K (MKL
 * Pardiso unless built with MKL_NO_OWN_PRODUCTS, cuDSS from
 * CUDSS_MIN_PRODUCT_NNZ: the rule each states up front) measures K'.
 *  (g) res is the residual with the K of the solver's rule; the two K's
 *      residuals differ by at least 10 LSS_TOL ||b|| (asserted, so the rule
 *      is observable);
 *  (h) direct solvers with SCS's products: the correction form from the
 *      plain solution reduces the residual at least twofold (the warm start
 *      was used; a step contracts by about LSS_DELTA / (1 + LSS_DELTA)); a
 *      warm start 1e6 ||x|| away gives the plain residual (the fallback);
 *      refinement with at most 1, 2 and 5 steps gives strictly decreasing
 *      residuals, at most half and at most 1e-2 of the plain one (steps are
 *      taken, kept only when better, and max_refine_steps is honoured); tol
 *      twice the plain residual takes no step, half of it exactly one;
 *      max_refine_steps 0 takes none;
 *  (i) CG under the mismatch: restarts solve with K' again, so they cannot
 *      remove it; the residual with 5 restarts is never above that with none;
 *  (j) solvers with their own K (no mismatch to refine): tol twice the plain
 *      residual, and refinement with tol 0, still give accurate solutions
 *      whose res matches, through the solver's own multiplies.
 */

/* LSS_TOL: the bar of every accuracy and residual check, relative to
 * ||b||_inf (see above). LSS_DELTA: the relative change of R in phase 2,
 * far above the rounding of a residual and small enough that refinement
 * contracts fast. LSS_FAR: how far (times ||x||) the far warm start of the
 * direct solvers is. LSS_CG_TOL: the CG solvers' tolerance. LSS_CG_FAR: how
 * far the warm start of the CG restart case is: eps ||K s|| must be well
 * above 100 LSS_CG_TOL. */
#define LSS_TOL (1e-3)
#define LSS_DELTA (1e-1)
#define LSS_FAR (1e6)
#ifdef SFLOAT
#define LSS_CG_TOL (1e-4)
#define LSS_CG_FAR (1e5)
#else
#define LSS_CG_TOL (1e-10)
#define LSS_CG_FAR (1e10)
#endif

static void lss_free_matrix(ScsMatrix *M) {
  if (M) {
    scs_free(M->p);
    scs_free(M->i);
    scs_free(M->x);
    scs_free(M);
  }
}

/* A (m x n, two thirds of the entries, uniform in [-1, 1] / 16), P upper
 * triangular tridiagonal and diagonally dominant (PSD), diag_r in
 * [0.5, 1.5] */
static void lss_data(scs_int n, scs_int m, ScsMatrix **A_out,
                     ScsMatrix **P_out, scs_float *diag_r) {
  scs_int i, j, nz = 0;
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsMatrix *P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  A->m = m;
  A->n = n;
  A->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  A->i = (scs_int *)scs_calloc(n * m, sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(n * m, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      if ((i + 3 * j) % 3 != 0 || i == j) {
        A->i[nz] = i;
        A->x[nz] = rand_scs_float() / 16.;
        nz++;
      }
    }
    A->p[j + 1] = nz;
  }
  P->m = n;
  P->n = n;
  P->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  P->i = (scs_int *)scs_calloc(2 * n, sizeof(scs_int));
  P->x = (scs_float *)scs_calloc(2 * n, sizeof(scs_float));
  nz = 0;
  for (j = 0; j < n; ++j) {
    if (j > 0) {
      P->i[nz] = j - 1;
      P->x[nz] = 0.5 * rand_scs_float();
      nz++;
    }
    P->i[nz] = j;
    P->x[nz] = 2. + rand_scs_float();
    nz++;
    P->p[j + 1] = nz;
  }
  for (i = 0; i < n + m; ++i) {
    diag_r[i] = 1. + 0.5 * rand_scs_float();
  }
  *A_out = A;
  *P_out = P;
}

/* nonzeros of K's stored triangle, as SCS(form_kkt) builds it: the diagonal,
 * P's strictly upper entries and A */
static scs_int lss_kkt_nnz(const ScsMatrix *A, const ScsMatrix *P) {
  scs_int j, k, nnz = A->n + A->m + A->p[A->n];
  for (j = 0; j < P->n; ++j) {
    for (k = P->p[j]; k < P->p[j + 1]; ++k) {
      nnz += P->i[k] < j;
    }
  }
  return nnz;
}

/* whether this solver measures its residuals with its own copy of K, by the
 * rule each states up front (compile-time knowledge the tests share): MKL
 * Pardiso unless built with MKL_NO_OWN_PRODUCTS (its other exception, a K
 * beyond MKL's 32-bit integers in a DLONG build without BLAS64, is far larger
 * than any test), cuDSS from CUDSS_MIN_PRODUCT_NNZ nonzeros of the stored
 * triangle; the others with SCS's products */
static scs_int lss_own_products(scs_int kkt_nnz) {
  const char *method = scs_get_lin_sys_method();
  if (strstr(method, "mkl")) {
#ifdef MKL_NO_OWN_PRODUCTS
    return 0;
#else
    return 1;
#endif
  }
  if (strstr(method, "cuDSS")) {
    return kkt_nnz >= CUDSS_MIN_PRODUCT_NNZ;
  }
  return 0;
}

/* ||b - K x||_inf with diag_r, from the problem data (the test's own
 * reference, independent of the linear solver's code) */
static scs_float lss_ref_residual(const ScsMatrix *A, const ScsMatrix *P,
                                  const scs_float *diag_r, const scs_float *b,
                                  const scs_float *x, scs_float *r) {
  scs_int i, n = A->n, m = A->m;
  memset(r, 0, (n + m) * sizeof(scs_float));
  SCS(accum_by_atrans)(A, &(x[n]), r);
  SCS(accum_by_p)(P, x, r);
  SCS(accum_by_a)(A, x, &(r[n]));
  for (i = 0; i < n; ++i) {
    r[i] = b[i] - (r[i] + diag_r[i] * x[i]);
  }
  for (i = n; i < n + m; ++i) {
    r[i] = b[i] - (r[i] - diag_r[i] * x[i]);
  }
  return SCS(norm_inf)(r, n + m);
}

/* the workspace and vectors of one size */
typedef struct {
  ScsLinSysWork *p;
  const ScsMatrix *A, *P;
  const scs_float *dr;  /* R given to scs_init_lin_sys_work (SCS's K) */
  const scs_float *dr2; /* R the solver solves with (= dr in phase 1) */
  scs_int nm, own, cg;
  scs_float nb; /* ||b||_inf */
  scs_float *b, *x, *xp, *s, *r;
} LssCase;

/* x = the solve of c->b by scs_solve_lin_sys(s, tol, max_steps, &res);
 * returns its status */
static scs_int lss_solve(LssCase *c, scs_float *x, const scs_float *s,
                         scs_float tol, scs_int max_steps, scs_float *res) {
  memcpy(x, c->b, c->nm * sizeof(scs_float));
  return scs_solve_lin_sys(c->p, x, s, tol, max_steps, res);
}

/* res must be the residual of x with the K of the solver's rule, and x
 * accurate if accurate; prints both K's residuals. Returns an error message
 * or SCS_NULL. */
static const char *lss_check(LssCase *c, const char *what, const scs_float *x,
                             scs_float res, scs_int accurate) {
  scs_float ref_k, ref_k2, ref;
  ref_k = lss_ref_residual(c->A, c->P, c->dr, c->b, x, c->r);
  ref_k2 = lss_ref_residual(c->A, c->P, c->dr2, c->b, x, c->r);
  ref = c->own ? ref_k2 : ref_k;
  scs_printf("test_lin_sys_solve: %s: res %.3e, ||b - K x|| %.3e, ||b - K' "
             "x|| %.3e (||b|| %.3e)\n",
             what, (double)res, (double)ref_k, (double)ref_k2, (double)c->nb);
  mu_assert("test_lin_sys_solve: res is not the residual of the solver's K",
            ABS(res - ref) <= LSS_TOL * c->nb);
  if (accurate) {
    mu_assert("test_lin_sys_solve: the solve is not accurate",
              ref <= LSS_TOL * c->nb);
  }
  return SCS_NULL;
}

#define LSS_CALL(st, msg)                                                      \
  mu_assert("test_lin_sys_solve: " msg " failed", (st) == 0)

#define LSS_CHECK(call)                                                        \
  do {                                                                         \
    const char *lss_msg_ = (call);                                             \
    if (lss_msg_) {                                                            \
      return lss_msg_;                                                         \
    }                                                                          \
  } while (0)

/* phase 1, K as given; returns an error message or SCS_NULL */
static const char *lss_phase1(LssCase *c) {
  scs_int i, nm = c->nm;
  scs_float res, res_p, res0, xn;
  /* (a) plain */
  LSS_CALL(lss_solve(c, c->xp, SCS_NULL, LSS_CG_TOL, 0, &res_p),
           "plain solve");
  LSS_CHECK(lss_check(c, "(a) plain", c->xp, res_p, 1));
  xn = SCS(norm_inf)(c->xp, nm);
  /* (b) a warm start near the solution */
  for (i = 0; i < nm; ++i) {
    c->s[i] = c->xp[i] + 1e-2 * xn * rand_scs_float();
  }
  LSS_CALL(lss_solve(c, c->x, c->s, LSS_CG_TOL, 0, &res),
           "solve from a near warm start");
  LSS_CHECK(lss_check(c, "(b) near warm start", c->x, res, 1));
  if (!c->cg) {
    /* (c) a warm start far from the solution: the fallback */
    for (i = 0; i < nm; ++i) {
      c->s[i] = LSS_FAR * xn * rand_scs_float();
    }
    LSS_CALL(lss_solve(c, c->x, c->s, LSS_CG_TOL, 0, &res),
             "solve from a far warm start");
    LSS_CHECK(lss_check(c, "(c) far warm start", c->x, res, 1));
    mu_assert("test_lin_sys_solve: a direct solver used a far warm start "
              "instead of solving plainly",
              res <= 100. * res_p);
  } else {
    /* (d) CG restarts */
    for (i = 0; i < nm; ++i) {
      c->s[i] = LSS_CG_FAR * rand_scs_float();
    }
    LSS_CALL(lss_solve(c, c->x, c->s, LSS_CG_TOL, 0, &res0),
             "CG from a very far warm start");
    LSS_CHECK(lss_check(c, "(d) CG from a very far warm start", c->x, res0, 0));
    mu_assert("test_lin_sys_solve: premise: CG's tracked residual from the "
              "very far warm start did not drift from the true one",
              res0 > 100. * LSS_CG_TOL);
    LSS_CALL(lss_solve(c, c->x, c->s, LSS_CG_TOL, 2, &res),
             "CG with restarts from a very far warm start");
    LSS_CHECK(lss_check(c, "(d) the same, with 2 restarts", c->x, res, 0));
    mu_assert("test_lin_sys_solve: CG's restarts did not reach the tolerance",
              res <= 10. * LSS_CG_TOL && res < res0);
  }
  /* (e) refinement: never worse */
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, c->cg ? LSS_CG_TOL : 0., 3, &res),
           "refined solve");
  LSS_CHECK(lss_check(c, "(e) refined, at most 3 steps", c->x, res, 1));
  mu_assert("test_lin_sys_solve: refinement made the solve worse",
            res <= 2. * res_p);
  return SCS_NULL;
}

/* phase 2: the solver solves with dr2; returns an error message or
 * SCS_NULL */
static const char *lss_phase2(LssCase *c) {
  scs_int i, nm = c->nm;
  scs_float res, res1, res2, res_p, ref_k, ref_k2, xn;
  /* (g) plain: res with the K of the solver's rule */
  LSS_CALL(lss_solve(c, c->xp, SCS_NULL, LSS_CG_TOL, 0, &res_p),
           "plain solve after the update");
  LSS_CHECK(lss_check(c, "(g) plain after the update", c->xp, res_p, c->own));
  ref_k = lss_ref_residual(c->A, c->P, c->dr, c->b, c->xp, c->r);
  ref_k2 = lss_ref_residual(c->A, c->P, c->dr2, c->b, c->xp, c->r);
  mu_assert("test_lin_sys_solve: premise: the two K's residuals do not "
            "differ, so which K the solver measures with is not observable",
            ABS(ref_k - ref_k2) >= 10. * LSS_TOL * c->nb);
  xn = SCS(norm_inf)(c->xp, nm);
  if (c->own) {
    /* (j) the solver's own K follows the update: nothing to refine, but
     * the calls run through its multiplies */
    LSS_CALL(lss_solve(c, c->x, SCS_NULL, 2. * res_p, 5, &res),
             "own K: solve with a tol already met");
    LSS_CHECK(lss_check(c, "(j) tol twice the plain residual", c->x, res, 1));
    LSS_CALL(lss_solve(c, c->x, SCS_NULL, c->cg ? LSS_CG_TOL : 0., 3, &res),
             "own K: refined solve");
    LSS_CHECK(lss_check(c, "(j) refined, at most 3 steps", c->x, res, 1));
    scs_printf("test_lin_sys_solve: %s measures with its own K, which "
               "follows the update: no mismatch to refine\n",
               scs_get_lin_sys_method());
    return SCS_NULL;
  }
  if (c->cg) {
    /* (i) CG restarts solve with K' again: never worse, no gain */
    LSS_CALL(lss_solve(c, c->x, SCS_NULL, LSS_CG_TOL, 5, &res),
             "CG with restarts under the mismatch");
    LSS_CHECK(lss_check(c, "(i) CG, 5 restarts", c->x, res, 0));
    mu_assert("test_lin_sys_solve: CG's restarts made the solve worse",
              res <= res_p + LSS_TOL * c->nb);
    return SCS_NULL;
  }
  /* (h) direct solvers with SCS's products: K' solves, K residuals */
  LSS_CALL(lss_solve(c, c->x, c->xp, 0., 0, &res),
           "correction form from the plain solution");
  LSS_CHECK(lss_check(c, "(h) warm start = plain solution", c->x, res, 0));
  mu_assert("test_lin_sys_solve: the warm start near the solution was not "
            "used (the residual did not halve)",
            res <= 0.5 * res_p);
  for (i = 0; i < nm; ++i) {
    c->s[i] = LSS_FAR * xn * rand_scs_float();
  }
  LSS_CALL(lss_solve(c, c->x, c->s, 0., 0, &res), "far warm start");
  LSS_CHECK(lss_check(c, "(h) far warm start", c->x, res, 0));
  mu_assert("test_lin_sys_solve: a far warm start was not refused",
            res <= 2. * res_p);
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, 0., 1, &res1), "1 step");
  LSS_CHECK(lss_check(c, "(h) at most 1 step", c->x, res1, 0));
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, 0., 2, &res2), "2 steps");
  LSS_CHECK(lss_check(c, "(h) at most 2 steps", c->x, res2, 0));
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, 0., 5, &res), "5 steps");
  LSS_CHECK(lss_check(c, "(h) at most 5 steps", c->x, res, 0));
  mu_assert("test_lin_sys_solve: refinement steps did not contract, or more "
            "steps than allowed were taken",
            res_p > res1 && res1 > res2 && res1 <= 0.5 * res_p &&
                res <= 1e-2 * res_p);
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, 2. * res_p, 5, &res), "tol met");
  LSS_CHECK(lss_check(c, "(h) tol twice the plain residual", c->x, res, 0));
  mu_assert("test_lin_sys_solve: refined below a tol already met",
            res >= 0.5 * res_p);
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, 0.5 * res_p, 5, &res), "tol missed");
  LSS_CHECK(lss_check(c, "(h) tol half the plain residual", c->x, res, 0));
  mu_assert("test_lin_sys_solve: did not take exactly one step to a tol "
            "one step reaches",
            res <= 0.5 * res_p && res > 1e-2 * res_p);
  LSS_CALL(lss_solve(c, c->x, SCS_NULL, 0., 0, &res), "no steps");
  mu_assert("test_lin_sys_solve: refined with max_refine_steps 0",
            res >= 0.5 * res_p);
  return SCS_NULL;
}

/* both phases at one size; returns an error message or SCS_NULL */
static const char *lss_size(scs_int n, scs_int m) {
  scs_int i;
  ScsMatrix *A, *P;
  LssCase c;
  scs_float *diag_r = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  scs_float *diag_r2 = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  const char *fail = SCS_NULL;
  memset(&c, 0, sizeof(c));
  ran_start(17);
  lss_data(n, m, &A, &P, diag_r);
  c.nm = n + m;
  c.A = A;
  c.P = P;
  c.dr = c.dr2 = diag_r;
  c.own = lss_own_products(lss_kkt_nnz(A, P));
  c.cg = strstr(scs_get_lin_sys_method(), "indirect") != SCS_NULL;
  c.b = (scs_float *)scs_calloc(c.nm, sizeof(scs_float));
  c.x = (scs_float *)scs_calloc(c.nm, sizeof(scs_float));
  c.xp = (scs_float *)scs_calloc(c.nm, sizeof(scs_float));
  c.s = (scs_float *)scs_calloc(c.nm, sizeof(scs_float));
  c.r = (scs_float *)scs_calloc(c.nm, sizeof(scs_float));
  for (i = 0; i < c.nm; ++i) {
    c.b[i] = rand_scs_float();
  }
  c.nb = SCS(norm_inf)(c.b, c.nm);
  scs_printf("test_lin_sys_solve: %s, n %li, m %li, %li nonzeros in K, "
             "residuals with %s\n",
             scs_get_lin_sys_method(), (long)n, (long)m,
             (long)lss_kkt_nnz(A, P), c.own ? "its own K" : "SCS's products");
  c.p = scs_init_lin_sys_work(A, P, diag_r);
  if (!c.p) {
    fail = "test_lin_sys_solve: scs_init_lin_sys_work failed";
  } else {
    fail = lss_phase1(&c);
    if (!fail) {
      /* the solver solves with R' = R (1 + LSS_DELTA), in another array
       * (kept alive: CG keeps the pointer); SCS's products keep diag_r */
      for (i = 0; i < c.nm; ++i) {
        diag_r2[i] = diag_r[i] * (1. + LSS_DELTA);
      }
      c.dr2 = diag_r2;
      if (scs_update_lin_sys_diag_r(c.p, diag_r2) != 0) {
        fail = "test_lin_sys_solve: scs_update_lin_sys_diag_r failed";
      } else {
        fail = lss_phase2(&c);
      }
    }
  }
  scs_free_lin_sys_work(c.p);
  lss_free_matrix(A);
  lss_free_matrix(P);
  scs_free(diag_r);
  scs_free(diag_r2);
  scs_free(c.b);
  scs_free(c.x);
  scs_free(c.xp);
  scs_free(c.s);
  scs_free(c.r);
  return fail;
}

static const char *test_lin_sys_solve(void) {
  const char *fail = lss_size(300, 400); /* about 8e4 nonzeros in K */
  if (!fail) {
    fail = lss_size(100, 150); /* about 1e4 */
  }
  return fail;
}
