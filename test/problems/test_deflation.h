#include "eigcg.h"
#include "glbopts.h"
#include "linalg.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"
#include <string.h>

/*
 * eigCG deflation, shared by the indirect solvers (linsys/eigcg.h).
 *
 * test_eigcg_core drives the shared routine through a test-defined solver: a
 * dense K = Q diag(lambda) Q' of size TD_N with a random orthogonal Q, the
 * identity as preconditioner (so the filter's diagonal D is the identity as
 * well), host vectors and ops of its own, and a CG that calls the hooks where
 * the solvers do. The eigenvalues: TD_TAIL of them spread evenly in log scale
 * over [TD_LAMBDA_MIN, TD_BULK_MIN) (twice the previous each), the rest over
 * [TD_BULK_MIN, 1], so a condition number of 1e5 that deflating the tail
 * brings down to 1e2: the case deflation is for. A cold solve to 1e-12 runs
 * about 400 CG iterations, so the harvest goes through about ten thick
 * restarts, and keeps 24 to 30 vectors (seeds 1 to 6 and 42; with seed 42,
 * which the test uses, 24 with Apple's Accelerate and 29 with NVPL). (A
 * spectrum spread evenly in log scale over [1e-4, 1] instead kept 8 to 10
 * vectors and deflation saved 3 to 4% of the iterations; one with its tail
 * down to 1e-6 kept none: the filter is relative to theta.) In a build with
 * deflation (USE_LAPACK, double precision, DEFLATE_VECTORS > 0):
 *  (a) the harvest keeps at least one vector; each kept vector w passes the
 *      filter, which, with rho its Rayleigh quotient, implies ||K w - rho w||
 *      <= t / (1 - t) rho ||w|| for t = EIGCG_RITZ_RES_TOL (rho minimizes
 *      the residual over scalars and is within t theta of the Ritz value
 *      theta); and the smallest rho is within 10% of lambda_min (the next
 *      eigenvalue is twice it; a broken harvest gives Ritz values off by
 *      order one, or a filter that keeps none; measured: 4e-13 to 1.5e-12);
 *  (b) a cold solve of a new right-hand side to 1e-8 with the set takes at
 *      most 0.8 times the CG iterations of the same solve without it (a
 *      deflation that does nothing gives 1.0; measured: 0.27 to 0.28, 88
 *      iterations against 314 to 322),
 *      and its residual, recomputed, is at most 100 times the tolerance (a
 *      projection inconsistent with its images leaves CG's tracked residual
 *      off the true one by the size of the correction);
 *  (c) with images that do not match the set (AW scaled by 100) and a
 *      residual in the span of W, the projection grows the residual (about
 *      99 times: asserted above EIGCG_UNDO_GROWTH as the premise) and is
 *      undone: x and r come back to within 1e-8 of the size of the
 *      correction (rounding leaves about 1e-16 of it; a missing undo leaves
 *      all of it);
 *  (d) after SCS(eigcg_reset) the set is empty and the next cold solve
 *      harvests a new one, in the window the first harvest used, here
 *      filled with NaN first: a harvest reads only what it has written
 *      itself (linsys/eigcg.h), so it keeps a set again, every entry finite
 *      (a read of a stale entry gives NaN, which kills the harvest or the
 *      set); this is why allocating the window once at setup, rather than
 *      per harvest, leaves the CPU solver's iterates as they were;
 *  (e) ops whose alloc fails at its first, second, third or fourth call (W,
 *      AW, V, VS): SCS(eigcg_init) fails, leaves deflation off (max 0) and
 *      nothing allocated;
 *  (f) a plain solve that returns between SCS(eigcg_begin_solve) and
 *      SCS(eigcg_end_solve) (the zero right-hand side of the solvers), and
 *      (g) a CG that stops before SCS(eigcg_start) (initial residual below
 *      tol), leave no set and allocate nothing, and the next cold solve
 *      harvests.
 * In a build without deflation, init succeeds and leaves it off (max 0),
 * checked against the build flags, allocates nothing, and the hooks leave x
 * and r untouched.
 *
 * test_lin_sys_deflation exercises the linear solver's own ops (the GPU
 * solver's on the GPU) through scs_solve_lin_sys: a K whose reduced matrix is
 * ill-conditioned even after the Jacobi preconditioner (A = Sigma V' with V
 * random orthogonal and singular values spread in log scale over [1e-2, 1],
 * R_x = 1e-4, so a condition number of a few thousand, plus a small tridiagonal
 * P), a cold solve to CG_BEST_TOL (the harvest), 5 warm solves from perturbed
 * starts, 3 cold solves of new right-hand sides, then a change of R written
 * into the same array, as SCS does, and the same again. Every solution's
 * residual, recomputed from the data, must be at most TD_LSS_BAR. CG counts
 * are not public, so this checks that deflation keeps the solves correct,
 * not that it ran; test_eigcg_core checks the mechanism. Indirect solvers in
 * double precision only (deflation is off in single precision, and this
 * conditioning is beyond its accuracy).
 */

#define TD_N (300)
#define TD_TAIL (10)
#define TD_LAMBDA_MIN (1e-5)
#define TD_BULK_MIN (1e-2)
/* TD_LSS_BAR: the bar of test_lin_sys_deflation's residuals, relative to
 * ||b||_inf = 1. CG's tracked residual is at most tol (1e-8 or less here),
 * and the true one drifts from it by rounding, about eps ||K|| ||x|| times
 * the square root of the iterations (||x|| about 2e3 here: about 1e-11). A
 * projection whose images or coefficients are wrong, or a set kept across a
 * change of R, leaves an error of the size of the correction (about delta
 * R_x ||x||, about 2e-2, for the change of R below). */
#define TD_LSS_BAR (1e-6)
#define TD_LSS_DELTA (0.1)

/* the test's solver: K dense (n x n, by column), identity preconditioner */
typedef struct {
  scs_int n;
  scs_float *K;
  scs_int allocs;    /* calls to alloc */
  scs_int fail_from; /* alloc fails from this call on (-1: never) */
  scs_int live;      /* allocations not freed */
} TdCtx;

/* the allocations of the ops when deflation is on: W, AW, V and VS */
#define TD_ALLOCS (4)

static void td_matvec(const TdCtx *c, const scs_float *x, scs_float *y) {
  scs_int i, j, n = c->n;
  memset(y, 0, n * sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      y[i] += c->K[(size_t)j * n + i] * x[j];
    }
  }
}

static scs_float *td_alloc(ScsLinSysWork *w, scs_int rows, scs_int cols) {
  TdCtx *c = (TdCtx *)w;
  scs_float *v;
  scs_int call = c->allocs++;
  if (c->fail_from >= 0 && call >= c->fail_from) {
    return SCS_NULL;
  }
  v = (scs_float *)scs_calloc((size_t)rows * cols, sizeof(scs_float));
  if (v) {
    c->live++;
  }
  return v;
}

static void td_free(ScsLinSysWork *w, scs_float *v) {
  TdCtx *c = (TdCtx *)w;
  if (v) {
    c->live--;
    scs_free(v);
  }
}

static scs_int td_scale_copy(ScsLinSysWork *w, scs_float *dst,
                             const scs_float *src, scs_float a) {
  scs_int t;
  for (t = 0; t < ((TdCtx *)w)->n; ++t) {
    dst[t] = a * src[t];
  }
  return 0;
}

static scs_int td_copy(ScsLinSysWork *w, scs_float *dst, const scs_float *src,
                       scs_int cols) {
  memcpy(dst, src, (size_t)cols * ((TdCtx *)w)->n * sizeof(scs_float));
  return 0;
}

static scs_int td_mul_small(ScsLinSysWork *w, scs_float *C, const scs_float *A,
                            scs_int k, const scs_float *B, scs_int cols) {
  scs_int i, j, l, n = ((TdCtx *)w)->n;
  for (j = 0; j < cols; ++j) {
    scs_float *cj = C + (size_t)j * n;
    memset(cj, 0, n * sizeof(scs_float));
    for (l = 0; l < k; ++l) {
      const scs_float *al = A + (size_t)l * n;
      scs_float blj = B[j * k + l];
      for (i = 0; i < n; ++i) {
        cj[i] += al[i] * blj;
      }
    }
  }
  return 0;
}

static scs_int td_probe(ScsLinSysWork *w, const scs_float *v, scs_float *Kv,
                        scs_float theta, scs_float *rr, scs_float *mm) {
  TdCtx *c = (TdCtx *)w;
  scs_int t;
  scs_float r2 = 0., m2 = 0.;
  td_matvec(c, v, Kv);
  for (t = 0; t < c->n; ++t) {
    scs_float mw = theta * v[t]; /* D = I */
    r2 += (Kv[t] - mw) * (Kv[t] - mw);
    m2 += mw * mw;
  }
  *rr = r2;
  *mm = m2;
  return 0;
}

static scs_int td_gram(ScsLinSysWork *w, const scs_float *W,
                       const scs_float *AW, scs_int k, scs_float *G) {
  scs_int a, b, n = ((TdCtx *)w)->n;
  for (a = 0; a < k; ++a) {
    for (b = 0; b < k; ++b) {
      G[b * k + a] = SCS(dot)(W + (size_t)a * n, AW + (size_t)b * n, n);
    }
  }
  return 0;
}

static scs_int td_coeffs(ScsLinSysWork *w, const scs_float *W, scs_int k,
                         const scs_float *r, scs_float *c) {
  scs_int j, n = ((TdCtx *)w)->n;
  for (j = 0; j < k; ++j) {
    c[j] = SCS(dot)(W + (size_t)j * n, r, n);
  }
  return 0;
}

static scs_int td_project(ScsLinSysWork *w, const scs_float *W,
                          const scs_float *AW, scs_int k, const scs_float *c,
                          scs_float *x, scs_float *r, scs_int undo) {
  scs_int j, t, n = ((TdCtx *)w)->n;
  scs_float sg = undo ? -1. : 1.;
  for (j = 0; j < k; ++j) {
    for (t = 0; t < n; ++t) {
      x[t] += sg * c[j] * W[(size_t)j * n + t];
      r[t] -= sg * c[j] * AW[(size_t)j * n + t];
    }
  }
  return 0;
}

static scs_float td_norm(ScsLinSysWork *w, const scs_float *r) {
  return SCS(norm_inf)(r, ((TdCtx *)w)->n);
}

static const ScsEigcgOps td_ops = {td_alloc,     td_free,  td_scale_copy,
                                   td_copy,      td_mul_small, td_probe,
                                   td_gram,      td_coeffs, td_project,
                                   td_norm};

/* K = Q diag(lambda) Q' with Q orthogonal (Gram-Schmidt, twice, of a random
 * matrix) and lambda as above */
static scs_float *td_matrix(scs_int n) {
  scs_int i, j, k, pass;
  scs_float *Q = (scs_float *)scs_calloc((size_t)n * n, sizeof(scs_float));
  scs_float *K = (scs_float *)scs_calloc((size_t)n * n, sizeof(scs_float));
  scs_float *lam = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (i = 0; i < n * n; ++i) {
    Q[i] = rand_scs_float();
  }
  for (j = 0; j < n; ++j) {
    scs_float *qj = Q + (size_t)j * n, nrm;
    for (pass = 0; pass < 2; ++pass) {
      for (k = 0; k < j; ++k) {
        scs_float d = SCS(dot)(Q + (size_t)k * n, qj, n);
        SCS(add_scaled_array)(qj, Q + (size_t)k * n, n, -d);
      }
    }
    nrm = SCS(norm_2)(qj, n);
    SCS(scale_array)(qj, 1. / nrm, n);
    lam[j] = j < TD_TAIL
                 ? TD_LAMBDA_MIN * POWF(TD_BULK_MIN / TD_LAMBDA_MIN,
                                        (scs_float)j / TD_TAIL)
                 : TD_BULK_MIN * POWF(1. / TD_BULK_MIN, (scs_float)(j - TD_TAIL) /
                                                            (n - 1 - TD_TAIL));
  }
  for (j = 0; j < n; ++j) {
    for (i = j; i < n; ++i) {
      scs_float s = 0.;
      for (k = 0; k < n; ++k) {
        s += Q[(size_t)k * n + i] * lam[k] * Q[(size_t)k * n + j];
      }
      K[(size_t)j * n + i] = K[(size_t)i * n + j] = s;
    }
  }
  scs_free(Q);
  scs_free(lam);
  return K;
}

/* the test's CG on K x = b from zero to tol, with the hooks where the
 * solvers have them (without the projection if deflate is 0); returns the
 * iterations, and the products of the extraction in *probes */
static scs_int td_cg(TdCtx *c, ScsEigcg *e, const scs_float *b, scs_float *x,
                     scs_float tol, scs_int deflate, scs_int *probes) {
  scs_int i, its = 0, n = c->n;
  scs_float ztr, ztr_prev, alpha, beta;
  scs_float *r = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *p = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *Kp = (scs_float *)scs_calloc(n, sizeof(scs_float));
  SCS(eigcg_begin_solve)(e, 1);
  memset(x, 0, n * sizeof(scs_float));
  memcpy(r, b, n * sizeof(scs_float));
  if (deflate) {
    (void)SCS(eigcg_deflate)(e, x, r); /* the test's ops do not fail */
  }
  if (SCS(norm_inf)(r, n) >= tol) {
    /* z = r: the preconditioner is the identity */
    ztr = SCS(dot)(r, r, n);
    memcpy(p, r, n * sizeof(scs_float));
    SCS(eigcg_start)(e, r, ztr);
    for (i = 0; i < 10 * n; ++i) {
      td_matvec(c, p, Kp);
      alpha = ztr / SCS(dot)(p, Kp, n);
      SCS(add_scaled_array)(x, p, n, alpha);
      SCS(add_scaled_array)(r, Kp, n, -alpha);
      its = i + 1;
      if (SCS(norm_inf)(r, n) < tol) {
        break;
      }
      ztr_prev = ztr;
      ztr = SCS(dot)(r, r, n);
      if (ztr_prev == 0.) {
        break;
      }
      beta = ztr / ztr_prev;
      SCS(eigcg_step)(e, alpha, beta, r, ztr);
      SCS(scale_array)(p, beta, n);
      SCS(add_scaled_array)(p, r, n, 1.);
    }
  }
  *probes = SCS(eigcg_end_solve)(e);
  scs_free(r);
  scs_free(p);
  scs_free(Kp);
  return its;
}

/* ||b - K x||_inf */
static scs_float td_residual(const TdCtx *c, const scs_float *b,
                             const scs_float *x) {
  scs_int i;
  scs_float *r = (scs_float *)scs_calloc(c->n, sizeof(scs_float)), nrm;
  td_matvec(c, x, r);
  for (i = 0; i < c->n; ++i) {
    r[i] = b[i] - r[i];
  }
  nrm = SCS(norm_inf)(r, c->n);
  scs_free(r);
  return nrm;
}

/* whether this build deflates (glbopts.h: USE_LAPACK, not SFLOAT,
 * DEFLATE_VECTORS > 0), compile-time knowledge the test shares */
#if defined(USE_LAPACK) && !defined(SFLOAT) && DEFLATE_VECTORS > 0
#define TD_DEFLATES (1)
#else
#define TD_DEFLATES (0)
#endif

/* NaN in every window entry that SCS's code, or the test's ops, could read
 * stale: V, VS, T, E, ev and S (the outputs of LAPACK and BLAS calls, TS, H,
 * Gr and work, are written by those calls before anything reads them) */
static void td_fill_window_nan(ScsEigcg *e) {
  size_t i, win = (size_t)e->win, n = (size_t)e->n, nev2 = 2 * (size_t)e->max;
  const scs_float nan = NAN;
  for (i = 0; i < n * win; ++i) {
    e->V[i] = nan;
  }
  for (i = 0; i < n * nev2; ++i) {
    e->VS[i] = nan;
  }
  for (i = 0; i < win * win; ++i) {
    e->T[i] = nan;
    e->E[i] = nan;
  }
  for (i = 0; i < win; ++i) {
    e->ev[i] = nan;
  }
  for (i = 0; i < win * nev2; ++i) {
    e->S[i] = nan;
  }
}

/* entries of the set (W, AW and the Cholesky factor of W' AW) that are not
 * finite */
static scs_int td_set_nonfinite(const ScsEigcg *e) {
  size_t i, nk = (size_t)e->n * e->count;
  scs_int bad = 0;
  for (i = 0; i < nk; ++i) {
    bad += !(ABS(e->W[i]) < INFINITY) + !(ABS(e->AW[i]) < INFINITY);
  }
  for (i = 0; i < (size_t)e->count * e->count; ++i) {
    bad += !(ABS(e->G[i]) < INFINITY);
  }
  return bad;
}

/* the checks of a build with deflation; returns an error message or
 * SCS_NULL */
static const char *td_core_on(TdCtx *c, ScsEigcg *e, scs_float *b1,
                              scs_float *b2, scs_float *x, scs_float *r) {
  scs_int a, j, n = c->n, its0, its1, probes, kd, its_a, allocs, bad;
  scs_float rho_min = 0., res, scale, nx, nr, grow, err;
  const scs_float lmin = TD_LAMBDA_MIN;
  scs_float *Kw = (scs_float *)scs_calloc(n, sizeof(scs_float));
  const char *fail = SCS_NULL;
  /* (a) the harvest */
  its0 = td_cg(c, e, b1, x, 1e-12, 1, &probes);
  its_a = its0;
  kd = e->count;
  scs_printf("test_eigcg_core: (a) harvest: %li CG iterations, %li probed, "
             "%li kept of at most %li (window %li)\n",
             (long)its0, (long)probes, (long)kd, (long)e->max, (long)e->win);
  if (kd < 1) {
    scs_free(Kw);
    return "test_eigcg_core: (a) the harvest kept no vector";
  }
  for (a = 0; a < kd && !fail; ++a) {
    const scs_float *wa = e->W + (size_t)a * n;
    scs_float ww = SCS(dot)(wa, wa, n), rho, er;
    const scs_float t = EIGCG_RITZ_RES_TOL;
    td_matvec(c, wa, Kw);
    rho = SCS(dot)(wa, Kw, n) / ww;
    SCS(add_scaled_array)(Kw, wa, n, -rho);
    er = SCS(norm_2)(Kw, n) / (rho * SQRTF(ww));
    rho_min = a == 0 ? rho : MIN(rho_min, rho);
    if (!(er <= t / (1. - t))) {
      scs_printf("test_eigcg_core: (a) vector %li: Rayleigh quotient %.4e, "
                 "relative eigenresidual %.3e\n",
                 (long)a, (double)rho, (double)er);
      fail = "test_eigcg_core: (a) a kept vector fails the filter's bound";
    }
  }
  scs_free(Kw);
  if (fail) {
    return fail;
  }
  err = ABS(rho_min - lmin) / lmin;
  scs_printf("test_eigcg_core: (a) smallest Rayleigh quotient %.6e, "
             "lambda_min %.6e (relative error %.2e)\n",
             (double)rho_min, (double)lmin, (double)err);
  mu_assert("test_eigcg_core: (a) the set misses the smallest eigenvalue",
            err <= 0.1);
  /* (b) a new right-hand side, without and with the set */
  its0 = td_cg(c, e, b2, x, 1e-8, 0, &probes);
  its1 = td_cg(c, e, b2, x, 1e-8, 1, &probes);
  res = td_residual(c, b2, x);
  scs_printf("test_eigcg_core: (b) CG iterations to 1e-8: %li without the "
             "set, %li with it (ratio %.3f); residual %.3e\n",
             (long)its0, (long)its1, (double)its1 / its0, (double)res);
  mu_assert("test_eigcg_core: (b) deflation did not reduce the CG iterations",
            its1 <= 0.8 * its0);
  mu_assert("test_eigcg_core: (b) the deflated solve is not accurate",
            res <= 100. * 1e-8);
  mu_assert("test_eigcg_core: (b) a solve with a set harvested again",
            probes == 0 && e->count == kd);
  /* (c) images that do not match: AW scaled by 100, r = W(:, 0) */
  for (j = 0; j < kd * n; ++j) {
    e->AW[j] *= 100.;
  }
  memset(x, 0, n * sizeof(scs_float));
  memcpy(r, e->W, n * sizeof(scs_float));
  (void)SCS(eigcg_deflate)(e, x, r);
  /* from e->c, the coefficients it used: the premise (the residual the
   * projection gave, W(:, 0) - AW c, grew) and the size of the correction;
   * b2 is free now, as scratch */
  scale = 0.;
  memcpy(b2, e->W, n * sizeof(scs_float));
  for (j = 0; j < kd; ++j) {
    SCS(add_scaled_array)(b2, e->AW + (size_t)j * n, n, -e->c[j]);
    scale += ABS(e->c[j]) * (SCS(norm_inf)(e->W + (size_t)j * n, n) +
                             SCS(norm_inf)(e->AW + (size_t)j * n, n));
  }
  grow = SCS(norm_inf)(b2, n) / SCS(norm_inf)(e->W, n);
  nx = SCS(norm_inf)(x, n);
  SCS(add_scaled_array)(r, e->W, n, -1.);
  nr = SCS(norm_inf)(r, n);
  scs_printf("test_eigcg_core: (c) the projection would grow r %.1f times; "
             "after the undo ||x|| %.3e, ||r - r0|| %.3e (correction %.3e)\n",
             (double)grow, (double)nx, (double)nr, (double)scale);
  mu_assert("test_eigcg_core: (c) premise: the projection did not grow r",
            grow > EIGCG_UNDO_GROWTH);
  mu_assert("test_eigcg_core: (c) the projection was not undone",
            nx <= 1e-8 * scale && nr <= 1e-8 * scale);
  /* the true images back (the same products, so the same values) */
  for (j = 0; j < kd; ++j) {
    td_matvec(c, e->W + (size_t)j * n, e->AW + (size_t)j * n);
  }
  /* (d) reset, then a new harvest, in a window of NaN */
  SCS(eigcg_reset)(e);
  mu_assert("test_eigcg_core: (d) reset kept the set", e->count == 0);
  td_fill_window_nan(e);
  its1 = td_cg(c, e, b1, x, 1e-12, 1, &probes);
  bad = td_set_nonfinite(e);
  scs_printf("test_eigcg_core: (d) after the reset, in a window of NaN: %li "
             "CG iterations, %li probed, %li kept (the first harvest: %li, "
             "%li kept); %li entries of the set not finite\n",
             (long)its1, (long)probes, (long)e->count, (long)its_a, (long)kd,
             (long)bad);
  mu_assert("test_eigcg_core: (d) no new harvest after the reset",
            probes > 0 && e->count > 0);
  mu_assert("test_eigcg_core: (d) the harvest read a stale window entry",
            bad == 0);
  /* (f) a solve that returns between begin and end, then (g) a CG that
   * stops before it starts, then a harvest */
  allocs = c->allocs;
  SCS(eigcg_reset)(e);
  SCS(eigcg_begin_solve)(e, 1);
  memset(r, 0, n * sizeof(scs_float));
  its0 = td_cg(c, e, r, x, 1e-12, 1, &probes); /* b = 0: no iteration */
  mu_assert("test_eigcg_core: (f, g) a solve without iterations left a set",
            its0 == 0 && probes == 0 && e->count == 0);
  td_cg(c, e, b1, x, 1e-12, 1, &probes);
  mu_assert("test_eigcg_core: (f, g) no harvest after a solve without "
            "iterations",
            e->count > 0);
  mu_assert("test_eigcg_core: (f, g) memory allocated after setup",
            c->allocs == allocs && c->live == TD_ALLOCS);
  return SCS_NULL;
}

/* (e): alloc failures at init; returns an error message or SCS_NULL */
static const char *td_core_alloc_fail(TdCtx *c) {
  ScsEigcg e;
  scs_int k, status;
  for (k = 0; k < TD_ALLOCS; ++k) {
    c->allocs = 0;
    c->fail_from = k;
    status = SCS(eigcg_init)(&e, (ScsLinSysWork *)c, &td_ops, c->n);
    mu_assert("test_eigcg_core: (e) init succeeded without its memory",
              status != 0);
    mu_assert("test_eigcg_core: (e) a failed init left deflation on, memory "
              "allocated, or allocated after the failure",
              e.max == 0 && c->live == 0 && c->allocs == k + 1);
    SCS(eigcg_free)(&e); /* safe after a failed init */
    mu_assert("test_eigcg_core: (e) SCS(eigcg_free) after a failed init",
              c->live == 0 && e.max == 0);
  }
  c->fail_from = -1;
  return SCS_NULL;
}

static const char *test_eigcg_core(void) {
  TdCtx c;
  ScsEigcg e;
  scs_int i, n = TD_N, probes;
  scs_float *b1 = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *b2 = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *r = (scs_float *)scs_calloc(n, sizeof(scs_float));
  const char *fail = SCS_NULL;
  memset(&c, 0, sizeof(c));
  c.n = n;
  c.fail_from = -1;
  ran_start(42);
  for (i = 0; i < n; ++i) {
    b1[i] = rand_scs_float();
    b2[i] = rand_scs_float();
  }
  if (SCS(eigcg_init)(&e, (ScsLinSysWork *)&c, &td_ops, n) != 0) {
    fail = "test_eigcg_core: SCS(eigcg_init) failed";
  }
  scs_printf("test_eigcg_core: deflation %s (max %li), expected %s\n",
             e.max > 0 ? "on" : "off", (long)e.max, TD_DEFLATES ? "on" : "off");
  if (fail) {
    /* reported below */
  } else if ((e.max > 0) != TD_DEFLATES) {
    fail = "test_eigcg_core: deflation is not on exactly when the build "
           "flags allow it";
  } else if (e.max > 0) {
    c.K = td_matrix(n);
    fail = td_core_on(&c, &e, b1, b2, x, r);
    SCS(eigcg_free)(&e);
    if (!fail && c.live != 0) {
      fail = "test_eigcg_core: memory left after SCS(eigcg_free)";
    }
    if (!fail) {
      fail = td_core_alloc_fail(&c);
    }
  } else {
    /* off: no memory, no harvest, and the hooks leave x and r as they are */
    SCS(eigcg_begin_solve)(&e, 1);
    memcpy(x, b1, n * sizeof(scs_float));
    memcpy(r, b2, n * sizeof(scs_float));
    (void)SCS(eigcg_deflate)(&e, x, r);
    SCS(eigcg_start)(&e, r, 1.);
    SCS(eigcg_step)(&e, 1., 1., r, 1.);
    probes = SCS(eigcg_end_solve)(&e);
    if (c.live != 0 || c.allocs != 0 || e.harvest || e.count || probes ||
        memcmp(x, b1, n * sizeof(scs_float)) ||
        memcmp(r, b2, n * sizeof(scs_float))) {
      fail = "test_eigcg_core: deflation that is off did something";
    }
    SCS(eigcg_free)(&e);
  }
  scs_free(c.K);
  scs_free(b1);
  scs_free(b2);
  scs_free(x);
  scs_free(r);
  return fail;
}

/* ||b - K x||_inf, K = [R_x + P, A'; A, -R_y] with diag_r, from the data */
static scs_float td_kkt_residual(const ScsMatrix *A, const ScsMatrix *P,
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

/* A = Sigma V' (n x n, dense, CSC), P tridiagonal (upper triangle), diag_r:
 * R_x = 1e-4, R_y = 1 */
static void td_lss_data(scs_int n, ScsMatrix **A_out, ScsMatrix **P_out,
                        scs_float *diag_r) {
  scs_int i, j, k, pass, nz = 0;
  scs_float *V = (scs_float *)scs_calloc((size_t)n * n, sizeof(scs_float));
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsMatrix *P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  for (i = 0; i < n * n; ++i) {
    V[i] = rand_scs_float();
  }
  for (j = 0; j < n; ++j) {
    scs_float *vj = V + (size_t)j * n;
    for (pass = 0; pass < 2; ++pass) {
      for (k = 0; k < j; ++k) {
        scs_float d = SCS(dot)(V + (size_t)k * n, vj, n);
        SCS(add_scaled_array)(vj, V + (size_t)k * n, n, -d);
      }
    }
    SCS(scale_array)(vj, 1. / SCS(norm_2)(vj, n), n);
  }
  /* A(i, j) = sigma_i V(j, i) */
  A->m = A->n = n;
  A->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  A->i = (scs_int *)scs_calloc((size_t)n * n, sizeof(scs_int));
  A->x = (scs_float *)scs_calloc((size_t)n * n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      scs_float sigma = POWF(1e-2, 1. - (scs_float)i / (n - 1));
      A->i[nz] = i;
      A->x[nz] = sigma * V[(size_t)i * n + j];
      nz++;
    }
    A->p[j + 1] = nz;
  }
  P->m = P->n = n;
  P->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  P->i = (scs_int *)scs_calloc(2 * n, sizeof(scs_int));
  P->x = (scs_float *)scs_calloc(2 * n, sizeof(scs_float));
  nz = 0;
  for (j = 0; j < n; ++j) {
    if (j > 0) {
      P->i[nz] = j - 1;
      P->x[nz] = 0.5e-4 * rand_scs_float();
      nz++;
    }
    P->i[nz] = j;
    P->x[nz] = 2e-4;
    nz++;
    P->p[j + 1] = nz;
  }
  for (i = 0; i < n; ++i) {
    diag_r[i] = 1e-4;
    diag_r[n + i] = 1.;
  }
  scs_free(V);
  *A_out = A;
  *P_out = P;
}

static void td_free_matrix(ScsMatrix *M) {
  if (M) {
    scs_free(M->p);
    scs_free(M->i);
    scs_free(M->x);
    scs_free(M);
  }
}

/* one round (the harvest, warm solves, cold solves); returns an error
 * message or SCS_NULL */
static const char *td_lss_round(ScsLinSysWork *p, const ScsMatrix *A,
                                 const ScsMatrix *P, const scs_float *diag_r,
                                 scs_int round, scs_float *b, scs_float *x,
                                 scs_float *s, scs_float *r) {
  scs_int i, k, nm = A->n + A->m;
  scs_float res, xn, worst = 0.;
  for (i = 0; i < nm; ++i) {
    b[i] = rand_scs_float();
  }
  /* the cold solve to CG_BEST_TOL: the harvest */
  memcpy(x, b, nm * sizeof(scs_float));
  mu_assert("test_lin_sys_deflation: solve failed",
            scs_solve_lin_sys(p, x, SCS_NULL, CG_BEST_TOL, 0, SCS_NULL) == 0);
  res = td_kkt_residual(A, P, diag_r, b, x, r);
  xn = SCS(norm_inf)(x, nm);
  scs_printf("test_lin_sys_deflation: round %li: cold solve to %.0e: "
             "residual %.3e, ||x|| %.3e\n",
             (long)round, (double)CG_BEST_TOL, (double)res, (double)xn);
  mu_assert("test_lin_sys_deflation: the harvesting solve is not accurate",
            res <= TD_LSS_BAR);
  /* warm solves from perturbed starts */
  for (k = 0; k < 5; ++k) {
    for (i = 0; i < nm; ++i) {
      s[i] = x[i] + 1e-2 * xn * rand_scs_float();
    }
    memcpy(x, b, nm * sizeof(scs_float));
    mu_assert("test_lin_sys_deflation: solve failed",
              scs_solve_lin_sys(p, x, s, 1e-8, 0, SCS_NULL) == 0);
    worst = MAX(worst, td_kkt_residual(A, P, diag_r, b, x, r));
  }
  scs_printf("test_lin_sys_deflation: round %li: 5 warm solves to 1e-8: "
             "worst residual %.3e\n",
             (long)round, (double)worst);
  mu_assert("test_lin_sys_deflation: a warm deflated solve is not accurate",
            worst <= TD_LSS_BAR);
  /* cold solves of new right-hand sides */
  worst = 0.;
  for (k = 0; k < 3; ++k) {
    for (i = 0; i < nm; ++i) {
      b[i] = rand_scs_float();
    }
    memcpy(x, b, nm * sizeof(scs_float));
    mu_assert("test_lin_sys_deflation: solve failed",
              scs_solve_lin_sys(p, x, SCS_NULL, 1e-10, 0, SCS_NULL) == 0);
    worst = MAX(worst, td_kkt_residual(A, P, diag_r, b, x, r));
  }
  scs_printf("test_lin_sys_deflation: round %li: 3 cold solves to 1e-10: "
             "worst residual %.3e\n",
             (long)round, (double)worst);
  mu_assert("test_lin_sys_deflation: a cold deflated solve is not accurate",
            worst <= TD_LSS_BAR);
  return SCS_NULL;
}

static const char *test_lin_sys_deflation(void) {
  scs_int i, n = 200, nm = 2 * n;
  ScsMatrix *A, *P;
  ScsLinSysWork *p;
  scs_float *diag_r, *b, *x, *s, *r;
  const char *fail = SCS_NULL;
#ifdef SFLOAT
  return mu_skipped; /* see above */
#endif
  if (!strstr(scs_get_lin_sys_method(), "indirect")) {
    return mu_skipped; /* no CG, no deflation */
  }
  diag_r = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  b = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  x = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  s = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  r = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  ran_start(7);
  td_lss_data(n, &A, &P, diag_r);
  scs_printf("test_lin_sys_deflation: %s, n %li, m %li\n",
             scs_get_lin_sys_method(), (long)n, (long)n);
  p = scs_init_lin_sys_work(A, P, diag_r);
  if (!p) {
    fail = "test_lin_sys_deflation: scs_init_lin_sys_work failed";
  } else {
    fail = td_lss_round(p, A, P, diag_r, 1, b, x, s, r);
    if (!fail) {
      /* a new R in the same array, as SCS writes it: the set must go */
      for (i = 0; i < nm; ++i) {
        diag_r[i] *= 1. + TD_LSS_DELTA;
      }
      if (scs_update_lin_sys_diag_r(p, diag_r) != 0) {
        fail = "test_lin_sys_deflation: scs_update_lin_sys_diag_r failed";
      } else {
        fail = td_lss_round(p, A, P, diag_r, 2, b, x, s, r);
      }
    }
  }
  scs_free_lin_sys_work(p);
  td_free_matrix(A);
  td_free_matrix(P);
  scs_free(diag_r);
  scs_free(b);
  scs_free(x);
  scs_free(s);
  scs_free(r);
  return fail;
}
