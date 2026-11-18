#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"
#include <string.h>

/*
 * The optional products of a linear solver with its own copy of the KKT
 * matrix K = [R_x + P, A'; A, -R_y] (scs_lin_sys_residual, scs_lin_sys_keep
 * and scs_lin_sys_correct in linsys.h), against SCS's own product from the
 * problem data. A solver that does not provide them must say so
 * (SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED); the MKL and cuDSS solvers must
 * provide them. A solver that does must, for random b, x, R and a P with
 * off-diagonal entries:
 *  - return ||b - K x||_inf equal to SCS's to rounding (relative 1e-10 of
 *    the size of the terms, 1e-4 in single precision; a missing or
 *    transposed block, a sign error or a stale R gives a difference of order
 *    one);
 *  - after scs_lin_sys_keep of a solve's solution x1, give ||b - K x1|| for
 *    x == SCS_NULL, and from it scs_lin_sys_correct must return x1 +
 *    K^{-1}(b - K x1), the solution of K x = b (residual at most 1e-8
 *    relative, 0.1 in single precision; a wrong correction is off by order
 *    one), and equal SCS's own correction form x1 + K^{-1}(b - K x1) with
 *    the residual vector from the problem data, entry by entry (relative
 *    1e-8, 0.1 in single precision): the inf-norms above do not see a
 *    residual vector that is permuted or wrong in entries below its largest;
 *    the solver's residual vector itself is not visible, so it is compared
 *    through the solve;
 *  - do the same after scs_update_lin_sys_diag_r changed R.
 */

#define LSP_N (300)
#define LSP_M (400)

/* how close the solver's residual must be to SCS's (relative to the size of
 * the terms) and how accurate the correction form must be (relative to
 * ||b||); a bug gives errors of order one */
#ifdef SFLOAT
/* a single-precision MKL Pardiso solve of this K leaves a relative residual
 * of up to about 2.4e-3; MKL is not tested in single precision in CI, so
 * leave a wide margin (a bug still gives errors of order one) */
#define LSP_PROD_TOL (1e-4)
#define LSP_SOLVE_TOL (1e-1)
#else
#define LSP_PROD_TOL (1e-10)
#define LSP_SOLVE_TOL (1e-8)
#endif

static void lsp_free_matrix(ScsMatrix *M) {
  if (M) {
    scs_free(M->p);
    scs_free(M->i);
    scs_free(M->x);
    scs_free(M);
  }
}

/* random A (two thirds of the entries: K has about 8e4 nonzeros, above
 * the cuDSS solver's CUDSS_MIN_PRODUCT_NNZ), P upper triangular
 * tridiagonal and diagonally dominant (PSD), diag_r in [0.5, 1.5] times
 * 10^-k for a mix of magnitudes */
static void lsp_data(ScsMatrix **A_out, ScsMatrix **P_out, scs_float *diag_r) {
  scs_int n = LSP_N, m = LSP_M, i, j, nz = 0;
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
        A->x[nz] = rand_scs_float();
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
    diag_r[i] = (1. + 0.5 * rand_scs_float()) * (i % 3 == 0 ? 1e-3 : 1.);
  }
  *A_out = A;
  *P_out = P;
}

/* ||b - K x||_inf from the problem data (SCS's own product) */
static scs_float lsp_ref_residual(const ScsMatrix *A, const ScsMatrix *P,
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

/* checks the products at diag_r; returns an error message or SCS_NULL */
static const char *lsp_check(ScsLinSysWork *p, const ScsMatrix *A,
                             const ScsMatrix *P, const scs_float *diag_r,
                             scs_float *b, scs_float *x, scs_float *x1,
                             scs_float *r) {
  scs_int i, nm = A->n + A->m;
  scs_float hook, ref, scale, res, diff;
  for (i = 0; i < nm; ++i) {
    b[i] = rand_scs_float();
    x[i] = rand_scs_float();
  }
  /* the size of the terms whose rounding errors the residuals carry */
  scale = SCS(norm_inf)(b, nm) +
          (SCS(norm_inf)(A->x, A->p[A->n]) + SCS(norm_inf)(P->x, P->p[P->n]) +
           SCS(norm_inf)(diag_r, nm)) *
              (scs_float)nm * SCS(norm_inf)(x, nm);
  hook = scs_lin_sys_residual(p, b, x);
  ref = lsp_ref_residual(A, P, diag_r, b, x, r);
  scs_printf("test_lin_sys_products: ||b - K x||: solver %.17g, SCS %.17g\n",
             (double)hook, (double)ref);
  mu_assert("test_lin_sys_products: the residual differs from SCS's",
            ABS(hook - ref) <= LSP_PROD_TOL * scale);
  /* x1 = K^{-1} b0 (any solve), kept; then the correction form for b */
  memcpy(x1, x, nm * sizeof(scs_float));
  mu_assert("test_lin_sys_products: solve failed",
            scs_solve_lin_sys(p, x1, SCS_NULL, 1e-12) == 0);
  scs_lin_sys_keep(p, x1);
  hook = scs_lin_sys_residual(p, b, SCS_NULL);
  ref = lsp_ref_residual(A, P, diag_r, b, x1, r);
  scs_printf("test_lin_sys_products: ||b - K x0||: solver %.17g, SCS "
             "%.17g\n",
             (double)hook, (double)ref);
  mu_assert("test_lin_sys_products: the residual at the kept x0 differs from "
            "SCS's",
            ABS(hook - ref) <= LSP_PROD_TOL * scale);
  mu_assert("test_lin_sys_products: correct failed",
            scs_lin_sys_correct(p, x) == 0);
  res = lsp_ref_residual(A, P, diag_r, b, x, r);
  scs_printf("test_lin_sys_products: correction form, ||b - K x|| / ||b|| "
             "= %.3e\n",
             (double)(res / SCS(norm_inf)(b, nm)));
  mu_assert("test_lin_sys_products: the correction form does not solve "
            "K x = b",
            res <= LSP_SOLVE_TOL * SCS(norm_inf)(b, nm));
  /* SCS's own correction form in r: x1 + K^{-1}(b - K x1), the residual
   * from the problem data */
  lsp_ref_residual(A, P, diag_r, b, x1, r);
  mu_assert("test_lin_sys_products: solve failed",
            scs_solve_lin_sys(p, r, SCS_NULL, 1e-12) == 0);
  SCS(add_scaled_array)(r, x1, nm, 1.);
  diff = 0.;
  for (i = 0; i < nm; ++i) {
    diff = MAX(diff, ABS(x[i] - r[i]));
  }
  scs_printf("test_lin_sys_products: correction form, max |solver - SCS| / "
             "||x|| = %.3e\n",
             (double)(diff / MAX(1., SCS(norm_inf)(r, nm))));
  mu_assert("test_lin_sys_products: the correction form differs from SCS's",
            diff <= LSP_SOLVE_TOL * MAX(1., SCS(norm_inf)(r, nm)));
  return SCS_NULL;
}

static const char *test_lin_sys_products(void) {
  scs_int n = LSP_N, m = LSP_M, i;
  ScsMatrix *A, *P;
  ScsLinSysWork *p;
  scs_float *diag_r = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  scs_float *b = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  scs_float *x = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  scs_float *x1 = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  scs_float *r = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  const char *fail = SCS_NULL;
  const char *method = scs_get_lin_sys_method();
  scs_int expected = strstr(method, "cuDSS") || strstr(method, "mkl");
  ran_start(17);
  lsp_data(&A, &P, diag_r);
  p = scs_init_lin_sys_work(A, P, diag_r);
  if (!p) {
    fail = "test_lin_sys_products: scs_init_lin_sys_work failed";
  } else if (scs_lin_sys_residual(p, b, x) ==
             SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED) {
    scs_printf("test_lin_sys_products: %s does not provide the products\n",
               method);
    if (expected) {
      fail = "test_lin_sys_products: this solver should provide them";
    }
  } else {
    fail = lsp_check(p, A, P, diag_r, b, x, x1, r);
    if (!fail) {
      /* a new R: the solver's K must follow it */
      for (i = 0; i < n + m; ++i) {
        diag_r[i] *= (i % 2 == 0) ? 10. : 0.1;
      }
      if (scs_update_lin_sys_diag_r(p, diag_r) != 0) {
        fail = "test_lin_sys_products: scs_update_lin_sys_diag_r failed";
      } else {
        fail = lsp_check(p, A, P, diag_r, b, x, x1, r);
      }
    }
  }
  scs_free_lin_sys_work(p);
  lsp_free_matrix(A);
  lsp_free_matrix(P);
  scs_free(diag_r);
  scs_free(b);
  scs_free(x);
  scs_free(x1);
  scs_free(r);
  return fail;
}
