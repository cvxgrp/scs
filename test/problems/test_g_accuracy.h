#include "glbopts.h"
#include "linsys.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _WIN32
#include <io.h>
#define G_ACC_DUP _dup
#define G_ACC_DUP2 _dup2
#define G_ACC_FILENO _fileno
#define G_ACC_CLOSE _close
#else
#include <unistd.h>
#define G_ACC_DUP dup
#define G_ACC_DUP2 dup2
#define G_ACC_FILENO fileno
#define G_ACC_CLOSE close
#endif

/*
 * SCS solves K g = [c; -b], K = [R_x + P, A'; A, -R_y], once per
 * factorization and uses g in every iteration until the next one, so it
 * solves for g to about machine precision: it refines g until a step no
 * longer halves the residual ||[c; -b] - K g||_inf, measured against the size
 * of the terms it is computed from, || |[c; -b]| + |K| |g| ||_inf, and stops
 * with SCS_FAILED, asking for a bug report, if that solve fails at the setup
 * of a solve (glbopts.h, "Accuracy of the solve for g"; at a scale update it
 * goes back to the previous R instead, test_g_rollback at the end of this
 * file, as it does when the factorization fails there, test_factor_rollback).
 * test_g_accuracy checks both, for every linear solver the test binary is
 * built with, on a badly scaled QP that is
 * not normalized (so that K keeps the bad scaling), with that measure
 * computed here (test_g_accuracy_iteration below checks that g is as
 * accurate as the iteration needs):
 * 1. the g of the first factorization (a one-iteration solve), where a plain
 *    solve's is far above machine precision (printed);
 * 2. the g after a solve with scale updates (the update path: g solved again
 *    with the new R);
 * 3. b with a NaN entry (scs_update): SCS_FAILED, with a message that says
 *    so (the data, not a bug);
 * 4. the linear solver updated to R' = G_ACC_BAD_R R in an array of the
 *    test's (scs_update_lin_sys_diag_r), so that it solves with a K' far from
 *    the K that SCS measures g with, as a broken solver would (its g's
 *    measure is printed in the message): SCS_FAILED,
 *    with the message asking to open an issue. Skipped where the linear
 *    solver measures with its own K, which follows the update (MKL Pardiso
 *    unless built with MKL_NO_OWN_PRODUCTS; cuDSS from CUDSS_MIN_PRODUCT_NNZ
 *    nonzeros: the rule each states up front);
 * 5. the linear solver updated to that R with a NaN entry, so that its g is
 *    not finite (or it returns an error): SCS_FAILED, with the message
 *    asking to open an issue. Skipped with 4, for MKL Pardiso (what
 *    Pardiso does with a NaN pivot is not checked here), and for a solver
 *    that refuses the update.
 * The messages are checked unless printing is compiled out (NO_PRINTING).
 */

/* how badly the test data is scaled: single precision cannot solve the
 * double-precision instance to any accuracy, refined or not */
#ifdef SFLOAT
#define G_ACC_ROW_SCALE (1e2)
#define G_ACC_COL_SCALE (1e1)
#else
#define G_ACC_ROW_SCALE (1e4)
#define G_ACC_COL_SCALE (1e3)
#endif

/* The limit on g's measure, in both tests. SCS refines g until a step no
 * longer halves its residual, which ends at 0.3-0.4 machine epsilons of the
 * terms (median over random SOCPs; at most 1.7 in about 30,000 solves) and at
 * 0.06-0.97 in these tests (QDLDL, CG, the dense solver, Accelerate, MKL
 * Pardiso, GPU CG and cuDSS, in both precisions, on macOS, Linux and
 * Windows). A refinement that stops at a looser target is above 2: with
 * G_TOL_EPS 4, CG's g measured 2.3 (single precision, first factorization)
 * and 3.9 (double, after scale updates) in test_g_accuracy, and 2.9-3.9 in
 * test_g_accuracy_iteration. A plain solve is far above it (printed). The
 * limit is not tied to G_TOL_EPS, so that a looser G_TOL_EPS fails it. */
#define G_ACC_LIMIT (2 * SCS_MACH_EPS)
/* the factor on R of the solver in check 4: K' then solves for a g whose
 * measure against K is close to 1 (printed in the message), the most a g can
 * have */
#define G_ACC_BAD_R (1e12)

/* the badly scaled QP: random data (gen_random_prob_data), every fifth row
 * of A (and of b) times G_ACC_ROW_SCALE, every seventh column (and entry of
 * c) times G_ACC_COL_SCALE, and a diagonal P */
static void g_acc_data(ScsData *d, ScsCone *k, ScsSolution *opt_sol) {
  scs_int n = 80, m = 240, col_nnz = 9, nnz = n * col_nnz, i, j;
  k->z = 24;
  k->l = m - k->z;
  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, 4321);
  for (j = 0; j < n; ++j) {
    for (i = d->A->p[j]; i < d->A->p[j + 1]; ++i) {
      if (d->A->i[i] % 5 == 0) {
        d->A->x[i] *= G_ACC_ROW_SCALE;
      }
      if (j % 7 == 0) {
        d->A->x[i] *= G_ACC_COL_SCALE;
      }
    }
    if (j % 7 == 0) {
      d->c[j] *= G_ACC_COL_SCALE;
    }
  }
  for (i = 0; i < m; ++i) {
    if (i % 5 == 0) {
      d->b[i] *= G_ACC_ROW_SCALE;
    }
  }
  d->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->P->m = d->P->n = n;
  d->P->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  d->P->i = (scs_int *)scs_calloc(n, sizeof(scs_int));
  d->P->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    d->P->p[j + 1] = j + 1;
    d->P->i[j] = j;
    d->P->x[j] = (j % 3 == 0) ? 1e2 : 1e-3;
  }
}

/* ||[c; -b] - K g||_inf over || |[c; -b]| + |K| |g| ||_inf, the size of the
 * terms the residual is computed from (the measure SCS refines g by), both
 * computed here, in double whatever scs_float is. P is diagonal or absent in
 * these instances. */
static double g_acc_rel_res(const ScsWork *w, const scs_float *g) {
  scs_int i, j, p, n = w->d->n, m = w->d->m;
  const ScsMatrix *A = w->d->A, *P = w->d->P;
  double *r = (double *)scs_calloc(n + m, sizeof(double));
  double *t = (double *)scs_calloc(n + m, sizeof(double));
  double res = 0., size = 0., v;
  for (i = 0; i < n; ++i) {
    v = (double)w->diag_r[i] * g[i];
    r[i] = (double)w->d->c[i] - v;
    t[i] = fabs((double)w->d->c[i]) + fabs(v);
  }
  for (i = 0; i < m; ++i) {
    v = (double)w->diag_r[n + i] * g[n + i];
    r[n + i] = -(double)w->d->b[i] + v;
    t[n + i] = fabs((double)w->d->b[i]) + fabs(v);
  }
  for (j = 0; j < n; ++j) {
    for (p = P ? P->p[j] : 0; P && p < P->p[j + 1]; ++p) { /* diagonal */
      v = (double)P->x[p] * g[j];
      r[j] -= v;
      t[j] += fabs(v);
    }
    for (p = A->p[j]; p < A->p[j + 1]; ++p) {
      v = (double)A->x[p] * g[n + A->i[p]];
      r[j] -= v;
      t[j] += fabs(v);
      v = (double)A->x[p] * g[j];
      r[n + A->i[p]] -= v;
      t[n + A->i[p]] += fabs(v);
    }
  }
  for (i = 0; i < n + m; ++i) {
    res = MAX(res, fabs(r[i]));
    size = MAX(size, t[i]);
  }
  scs_free(r);
  scs_free(t);
  return res / size;
}

/* the same of a plain solve of the same system, to show what the refinement
 * had to fix (printed only: it depends on the linear solver) */
static double g_acc_unrefined(ScsWork *w) {
  scs_int i, n = w->d->n, m = w->d->m;
  scs_float *g0 = (scs_float *)scs_calloc(n + m, sizeof(scs_float));
  double res;
  memcpy(g0, w->d->c, n * sizeof(scs_float));
  for (i = 0; i < m; ++i) {
    g0[n + i] = -w->d->b[i];
  }
  scs_solve_lin_sys(w->p, g0, SCS_NULL, CG_BEST_TOL, 0, SCS_NULL);
  res = g_acc_rel_res(w, g0);
  scs_free(g0);
  return res;
}

/* scs_solve with its output captured; *found: whether the output contains
 * want and, unless it is SCS_NULL, want2, and does not contain unwanted
 * unless that is SCS_NULL (1 when printing is compiled out, or if the
 * capture cannot be set up: then the output goes to stdout as usual and only
 * the status is checked). Only the first G_ACC_CAPTURE_LEN - 1 bytes of the
 * output are read. */
#define G_ACC_CAPTURE_LEN (16384)
static scs_int g_acc_solve_capture(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                                   const char *want, const char *want2,
                                   const char *unwanted, scs_int *found) {
  scs_int status;
#if NO_PRINTING > 0
  *found = 1;
  return scs_solve(w, sol, info, 0);
#else
  FILE *tmp = tmpfile();
  int saved = -1;
  char buf[G_ACC_CAPTURE_LEN];
  size_t len;
  *found = 1;
  fflush(stdout);
  if (tmp) {
    saved = G_ACC_DUP(G_ACC_FILENO(stdout));
  }
  if (saved < 0 || G_ACC_DUP2(G_ACC_FILENO(tmp), G_ACC_FILENO(stdout)) < 0) {
    if (saved >= 0) {
      G_ACC_CLOSE(saved);
    }
    if (tmp) {
      fclose(tmp);
    }
    scs_printf("test_g_accuracy: cannot capture the output, the message is "
               "not checked\n");
    return scs_solve(w, sol, info, 0);
  }
  status = scs_solve(w, sol, info, 0);
  fflush(stdout);
  G_ACC_DUP2(saved, G_ACC_FILENO(stdout));
  G_ACC_CLOSE(saved);
  rewind(tmp);
  len = fread(buf, 1, sizeof(buf) - 1, tmp);
  buf[len] = '\0';
  fclose(tmp);
  scs_printf("%s", buf); /* the captured output, for the log */
  *found = strstr(buf, want) != SCS_NULL &&
           (!want2 || strstr(buf, want2) != SCS_NULL) &&
           (!unwanted || strstr(buf, unwanted) == SCS_NULL);
  return status;
#endif
}

static const char *test_g_accuracy(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs_float *b = SCS_NULL, *dr = SCS_NULL;
  ScsInfo info = {0};
  ScsWork *w;
  scs_int i, status, found, nnz;
  double res, res0;
  const char *fail = SCS_NULL, *method = scs_get_lin_sys_method();

  g_acc_data(d, k, opt_sol);
  scs_set_default_settings(stgs);
  stgs->normalize = 0;
  stgs->max_iters = 1;
  stgs->verbose = 0;

  w = scs_init(d, k, stgs);
  mu_assert("test_g_accuracy: scs_init returned NULL", w != SCS_NULL);
  /* 1. a solve computes g at its start (update_work); one iteration makes
   * no scale update, so g is the one of the first factorization */
  scs_solve(w, sol, &info, 0);
  res = g_acc_rel_res(w, w->g);
  res0 = g_acc_unrefined(w);
  scs_printf("test_g_accuracy: %s, first factorization: g relative residual "
             "%.2e (plain solve %.2e, limit %.2e)\n",
             method, res, res0, (double)G_ACC_LIMIT);
  if (info.scale_updates != 0) {
    fail = "test_g_accuracy: unexpected scale update in one iteration";
  } else if (!(res <= G_ACC_LIMIT)) {
    fail = "test_g_accuracy: g not accurate for the first factorization";
  }

  if (!fail) {
    /* 2. scale updates refactorize and solve for g again with the new R */
    w->stgs->max_iters = 2000;
    scs_solve(w, sol, &info, 0);
    res = g_acc_rel_res(w, w->g);
    res0 = g_acc_unrefined(w);
    scs_printf("test_g_accuracy: after %li iterations and %li scale updates, "
               "g relative residual %.2e (plain solve %.2e)\n",
               (long)info.iter, (long)info.scale_updates, res, res0);
    if (info.scale_updates < 1) {
      fail = "test_g_accuracy: no scale update, the second check is vacuous";
    } else if (!(res <= G_ACC_LIMIT)) {
      fail = "test_g_accuracy: g not accurate after scale updates";
    }
  }

  if (!fail) {
    /* 3. NaN in b: no g can be finite, and the message blames the data */
    b = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
    memcpy(b, d->b, d->m * sizeof(scs_float));
    b[1] = NAN;
    scs_update(w, b, SCS_NULL);
    w->stgs->max_iters = 10;
    status = g_acc_solve_capture(w, sol, &info, "b or c contains a NaN",
                                 SCS_NULL, SCS_NULL, &found);
    scs_printf("test_g_accuracy: NaN in b: status %li\n", (long)status);
    if (status != SCS_FAILED || info.status_val != SCS_FAILED) {
      fail = "test_g_accuracy: NaN in b did not fail";
    } else if (!found) {
      fail = "test_g_accuracy: NaN in b: no message about the data";
    }
    b[1] = d->b[1];
    scs_update(w, b, SCS_NULL);
  }

  nnz = d->n + d->m + d->A->p[d->n] + d->P->p[d->n]; /* K's stored triangle */
  if (!fail && (
#if MKL_NO_OWN_PRODUCTS <= 0
                   strstr(method, "mkl") ||
#endif
                   (strstr(method, "cuDSS") && nnz >= CUDSS_MIN_PRODUCT_NNZ))) {
    scs_printf("test_g_accuracy: %s measures with its own K, which follows "
               "the update: the broken-solver check is skipped\n",
               method);
  } else if (!fail) {
    /* 4. a solver that solves with K' = K(G_ACC_BAD_R R): the array stays
     * alive until scs_finish (CG keeps the pointer) */
    dr = (scs_float *)scs_calloc(d->n + d->m, sizeof(scs_float));
    for (i = 0; i < d->n + d->m; ++i) {
      dr[i] = w->diag_r[i] * G_ACC_BAD_R;
    }
    if (scs_update_lin_sys_diag_r(w->p, dr) != 0) {
      fail = "test_g_accuracy: update to the test's R failed";
    } else {
      status = g_acc_solve_capture(
          w, sol, &info, "it indicates a bug in SCS",
          "open an issue at https://github.com/cvxgrp/scs/issues", SCS_NULL,
          &found);
      scs_printf("test_g_accuracy: solver on %.0e R: status %li\n",
                 G_ACC_BAD_R, (long)status);
      if (status != SCS_FAILED || info.status_val != SCS_FAILED) {
        fail = "test_g_accuracy: a broken solve for g did not fail";
      } else if (!found) {
        fail = "test_g_accuracy: a broken solve for g: no message saying it "
               "is a bug and asking to open an issue";
      }
    }
    if (!fail && strstr(method, "mkl")) {
      scs_printf("test_g_accuracy: %s: the check of a g that is not finite "
                 "is skipped\n",
                 method);
    } else if (!fail) {
      /* 5. NaN in the solver's R: g0 is not finite (or the solver returns an
       * error), and either way SCS stops asking for a bug report */
      dr[0] = NAN;
      if (scs_update_lin_sys_diag_r(w->p, dr) != 0) {
        scs_printf("test_g_accuracy: %s refuses an R with a NaN entry: the "
                   "check of a g that is not finite is skipped\n",
                   method);
      } else {
        status = g_acc_solve_capture(
            w, sol, &info, "it indicates a bug in SCS",
            "open an issue at https://github.com/cvxgrp/scs/issues", SCS_NULL,
            &found);
        scs_printf("test_g_accuracy: solver on R with a NaN entry: status "
                   "%li\n",
                   (long)status);
        if (status != SCS_FAILED || info.status_val != SCS_FAILED) {
          fail = "test_g_accuracy: a g that is not finite did not fail";
        } else if (!found) {
          fail = "test_g_accuracy: a g that is not finite: no message saying "
                 "it is a bug and asking to open an issue";
        }
      }
    }
  }

  scs_finish(w);
  scs_free(b);
  scs_free(dr);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  return fail;
}

/*
 * test_g_accuracy_iteration: g as accurate as the iteration needs. An
 * iteration that uses a computed g is exact for the data K g, that is for
 * c - r_x and b + r_y (r = [c; -b] - K g), so g's residual is a floor under
 * SCS's residuals that no number of iterations removes. Relative to SCS's
 * tolerance it is largest where the precision is coarse and the rows are
 * dense (many terms per row, so a large size of the terms and of g's
 * rounding): in single precision at eps 1e-6, about eight machine epsilons,
 * on this random SOCP (100 variables, 130 rows, every entry of A nonzero, 50
 * linear rows and one second-order cone of 80). With G_TOL_EPS 4, the
 * refinement of single-precision CG stopped at 2.9-3.9 epsilons of g's
 * terms, and this problem ran to max_iters at every start tried (10,000
 * iterations, the dual residual stuck at 1e-5 to 3e-5); refined to the
 * halving stop (0.3-0.6 epsilons) it solved in 400 to 900 iterations at the
 * same 9 starts (initial scale moved by up to 4 ulps; Mac). The test runs in
 * every build: in double precision, and with the direct solvers in single,
 * g was accurate enough either way. Checked: g's measure at the end at most
 * G_ACC_LIMIT, in every build; in double precision, also that the solve
 * finished within G_ITER_MAX_ITERS.
 *
 * In single precision whether the solve finished is printed, not checked,
 * because it does not tell an inaccurate g from the iteration's own
 * rounding. Built without FMA contraction (-march=x86-64 or -v2, or
 * -ffp-contract=off), single precision at eps 1e-6 can stay above the
 * stopping thresholds for thousands of iterations with g accurate (its
 * measure at most 1.1 epsilons), the dual residual at 2.8e-5 (median)
 * against 2.5e-6 to 5.7e-6 where solved runs stop: up to 19,225 iterations
 * over 21 starting scales each with QDLDL and CG on the Mac, and past 20,000
 * at some starts with every linear solver on x86-64 hosts (54,400 at the
 * test's own start with the dense solver and OpenBLAS's AVX-512 kernels).
 * Builds with FMA contraction took at most 1,000. With G_TOL_EPS 4, CG ran
 * to 20,000 at each of the 21 starts, with FMA contraction or without, so
 * without it the cap does not separate the two. g's measure does, in every
 * build: with G_TOL_EPS 4, CG's g measured 2.8 to 4.1 epsilons at those
 * starts, against at most 1.1 with g refined to the halving stop.
 */
#define G_ITER_N (100)
#define G_ITER_SEED (4)
/* eps_abs = eps_rel, set with scs_update_settings after scs_init, which in
 * single-precision builds of this suite raises every tolerance to at least
 * 1e-4 (problem_utils.h): at 1e-4 the floor is far below the tolerance */
#define G_ITER_EPS (1e-6)
/* Single precision at 1e-6 on this generator is path-chaotic even with g
 * refined to the halving stop: over seeds 1-40, the most iterations were
 * 6,050 with CG and 13,650 with QDLDL (one seed each above 5,000; seed 4
 * took 400-900 with CG and 225-500 with QDLDL at 9 starts). 20,000, the
 * single-precision cap of the A/B, leaves room for a runner whose rounding
 * takes another path. The cap costs nothing while the test passes; with the
 * floor above the tolerance no number of iterations would do (G_TOL_EPS 4,
 * single-precision CG: still at max_iters after 20,000, g at 3.4 epsilons of
 * its terms, about a minute on the Mac). In single precision, where the test
 * does not check that the solve finished (see above), the cap only bounds
 * the time: a build without FMA contraction can run all of it, about 50
 * seconds with CG on the Mac. */
#define G_ITER_MAX_ITERS (20000)
/* whether the test checks that the solve finished within G_ITER_MAX_ITERS:
 * not in single precision (see above) */
#ifdef SFLOAT
#define G_ITER_CHECK_SOLVED (0)
#else
#define G_ITER_CHECK_SOLVED (1)
#endif

static const char *test_g_accuracy_iteration(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  ScsWork *w;
  scs_int exitflag, n = G_ITER_N, m = 13 * G_ITER_N / 10;
  double res;
  const char *fail = SCS_NULL;

  d->n = n;
  d->m = m;
  k->qsize = 1;
  k->q = (scs_int *)scs_calloc(1, sizeof(scs_int));
  k->q[0] = 4 * n / 5;
  k->l = m - k->q[0];
  /* col_nnz = m: every entry of A is nonzero */
  gen_random_prob_data(m * n, m, d, k, opt_sol, G_ITER_SEED);
  scs_set_default_settings(stgs);
  stgs->verbose = 0;
  w = scs_init(d, k, stgs);
  stgs->eps_abs = stgs->eps_rel = G_ITER_EPS;
  stgs->max_iters = G_ITER_MAX_ITERS;
  if (w == SCS_NULL) {
    fail = "test_g_accuracy_iteration: scs_init returned NULL";
  } else if (scs_update_settings(w, stgs) != 0) {
    fail = "test_g_accuracy_iteration: scs_update_settings failed";
  } else {
    exitflag = scs_solve(w, sol, &info, 0);
    res = g_acc_rel_res(w, w->g);
    scs_printf("test_g_accuracy_iteration: %s, eps %.0e: exit %li after %li "
               "iterations (%li scale updates), res_pri %.2e, res_dual "
               "%.2e; last g relative residual %.2e (limit %.2e)\n",
               scs_get_lin_sys_method(), G_ITER_EPS, (long)exitflag,
               (long)info.iter, (long)info.scale_updates, info.res_pri,
               info.res_dual, res, (double)G_ACC_LIMIT);
    if (G_ITER_CHECK_SOLVED && exitflag != SCS_SOLVED) {
      fail = "test_g_accuracy_iteration: not solved within max_iters "
             "(g not as accurate as the iteration needs?)";
    } else if (!(res <= G_ACC_LIMIT)) {
      fail = "test_g_accuracy_iteration: g not accurate";
    }
  }
  scs_finish(w);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  return fail;
}

/*
 * test_g_rollback: a scale update whose solve for g fails is rolled back
 * (scale_update_rollback in src/scs.c): SCS goes back to the previous R and
 * continues instead of stopping with SCS_FAILED, which only the setup of a
 * solve still does (test_g_accuracy, checks 4 and 5).
 *
 * How the solve for g is made to fail at a scale update and nowhere else.
 * Checks 4 and 5 above, and test_refine_perturbed, give the linear solver an
 * array of the test's, against linsys.h's rule that a new R is written into
 * the array given at init, so that it solves with a K' while SCS's products
 * keep measuring K. A scale update undoes that (SCS passes its own diag_r
 * again), which is why those tests turn adaptive_scale off. This test does
 * the mirror image: after scs_init it points w->diag_r at an array of its
 * own with the same values. SCS then writes every new R into the test's
 * array and hands that array to the linear solver, while SCS's products (the
 * residual of the shared routine, which the check of g uses) keep measuring
 * with the array given at init, which keeps the R of the setup. At the setup
 * the two agree, and the check passes. After an attempted update the solver
 * solves with the new R and the check measures with the old one, as if the
 * solver had failed at the new R. The rollback writes the old R back, and
 * the two agree again. No hook in the library.
 *
 * Why the scale setting is G_RB_SCALE. Any scalar update from a scale below
 * MIN_SCALE_VALUE (1e-6) is clamped to at least that, so it divides every
 * entry of R_y by at least 1e-6 / G_RB_SCALE. The residual the check measures
 * in the rows of y is then (R_old - R_new) g_y, that many times R_new g_y,
 * while the terms it is divided by are of the size of R_new g_y: the check
 * fails by a wide margin whatever the linear solver. The scale setting may be
 * any positive number (scale.rst: the limits apply to the updates only).
 *
 * What is checked. Every attempted update is rolled back, so R never
 * changes, and the solve must go on as if no update had been attempted, that
 * is as the same solve with adaptive_scale off, which never attempts one.
 * Not bit for bit: where the arithmetic is not reproducible from one solve to
 * the next, two solves differ by rounding, and this fixture amplifies it (MKL
 * in CI, 1000 iterations: pobj -624644852.7 against -624644853.6). With MKL
 * that is heap alignment, not threading: two solves in one process get their
 * arrays at other addresses (and the solve with adaptive_scale off allocates
 * no per-row multipliers), and MKL's BLAS and Pardiso, without conditional
 * numerical reproducibility, round according to the alignment. cuDSS and GPUs
 * are not reproducible either. Two runs (the instance of test_g_accuracy,
 * g_acc_data, with the default settings otherwise, verbose off, with the
 * output captured), each against the same solve with adaptive_scale off:
 * 1. G_RB_ITERS iterations with adaptive_diag_scale off, so that every
 *    attempt is a scalar one: all G_RB_ITERS iterations run (a failed
 *    restore stops the solve at once, with SCS_FAILED), at least one
 *    rollback, at most one per RESCALING_MIN_ITERS iterations (the pacing),
 *    no scale update kept, the scale still G_RB_SCALE, the rollback's
 *    WARNING line in the output (printed whatever verbose is) with the
 *    residual's ratio to the size of its terms, and no message asking for a
 *    bug report. Not the reference's status or iterates: after 1000
 *    iterations, rounding alone (as below) moved CG's iterates by up to 2 and
 *    changed the status of one of the two solves in 5 of 192 runs (double
 *    precision).
 * 2. G_RB_ITERS_DIAG iterations with adaptive_diag_scale on: the same, the
 *    per-row multipliers back at 1 (the attempt at iteration
 *    RESCALING_MIN_ITERS also installs them), the reference's status, and x,
 *    y and s within G_RB_TOL of the reference's: the largest over x, y and s
 *    of ||a - b||_inf / max(||a||_inf, ||b||_inf) (g_rb_dev).
 * Skipped where the linear solver measures with its own K, which follows the
 * update, as check 4 is.
 *
 * Not covered: a failed restore (the restored R is consistent by
 * construction here; test_factor_rollback below makes the restore's
 * factorization fail, and the failure of its solve for g is the setup's path,
 * checks 4 and 5); the cases where CG's deflation basis does not come back
 * bit for bit (a later solve on the same workspace, or a solve for g that
 * kept no vector); and the cone
 * cache: a rollback that also reset it (SCS(reset_cone_cache)) passes, since
 * the cache holds only the warm starts of the logdet cone's Newton solves
 * (USE_SPECTRAL_CONES builds) and this instance has no such cone. A PSD cone
 * would not help: without USE_SPECTRAL_CONES the reset does nothing.
 */
/* G_RB_SCALE: the scale setting (see above). G_RB_ITERS: the long run, room
 * for 9 attempts, one per RESCALING_MIN_ITERS. G_RB_ITERS_DIAG: one attempt,
 * at iteration RESCALING_MIN_ITERS, and the iterates compared 50 iterations
 * later. G_RB_TOL: the largest g_rb_dev of run 2.
 *
 * Measured on the Mac (QDLDL, CG and dense): an attempt at every check where
 * one is allowed (iterations 100 to 900), every one rolled back, the check's
 * ratio far above G_FAIL_REL_RES (0.5). In double precision G_RB_SCALE is
 * 1e-12: every attempt clamped from 1e-12 to 1e-6, ratio 5e5 to 6e6. In single
 * precision it is 1e-7: at 1e-12, K is far beyond single precision, and
 * rounding alone (as below) moved run 2's iterates by up to 2e-2 and changed
 * its status in 11 of 768 runs; at 1e-7 the attempts move the scale 46- to
 * 183-fold (to 4.6e-6 to 1.8e-5), ratio 22 to 88. Either way the iterates
 * diverge or stall: the runs are a fixture for the rollback, not a solve.
 *
 * G_RB_TOL sits between the two measured sizes of g_rb_dev in run 2:
 * - rounding: every BLAS dot and norm, or every entry of every linear
 *   solve's answer, or both, moved by -1, 0 or +1 ulp at random (64 seeds
 *   each, per solver and precision), OpenBLAS 0.3.34 threaded against
 *   Accelerate, one entry of b or c or the scale moved by one ulp: at most
 *   1.8e-10 in double and 2.1e-3 in single, and no change of the status,
 *   the rollbacks or the multipliers (nor of run 1's rollbacks);
 * - a broken rollback (the review's mutants): the acceleration reset at the
 *   rollback, 0.096 (CG) to 0.20 (QDLDL, dense) in double and 0.056 in
 *   single; the linear solver not given the previous R back, 1.3 (CG, double)
 *   and 0.076 (CG, single), SCS_FAILED with QDLDL and dense; g not solved
 *   again at the restore, SCS_FAILED or a certificate status in double, 0.96
 *   in single; the scale or the multipliers not put back, SCS_FAILED; no
 *   pacing after a rollback, 36 rollbacks in run 1.
 * So G_RB_TOL is 1e-3 in double, 96 times below the smallest effect and 5.6e6
 * times above rounding, and 1e-2 in single, where the two are only 27 times
 * apart: 4.8 times above rounding and 5.6 times below the smallest effect.
 * Where the arithmetic is reproducible from one solve to the next (on the
 * Mac: Accelerate's or OpenBLAS's BLAS at their default threading, QDLDL, CG,
 * dense and Accelerate's solver), the two runs are the same bits and
 * g_rb_dev is 0. */
#ifdef SFLOAT
#define G_RB_SCALE (1e-7)
#define G_RB_TOL (1e-2)
#else
#define G_RB_SCALE (1e-12)
#define G_RB_TOL (1e-3)
#endif
#define G_RB_ITERS (1000)
#define G_RB_ITERS_DIAG (150)

/* One solve of g_acc_data's instance at scale G_RB_SCALE and at most
 * max_iters iterations, adaptive_scale and adaptive_diag_scale as given, with
 * SCS's diag_r moved to an array of the test's as described above. Returns
 * the status; *rollbacks is w->scale_rollbacks, *mults_one whether every
 * per-row multiplier is 1 at the end (1 if there are none), *found whether
 * the output has the rollback's WARNING line, for a solve for g that failed
 * its check, and no message asking for a bug report
 * (see g_acc_solve_capture), or -1 for every output if setup failed. */
static scs_int g_rb_solve(scs_int adaptive, scs_int diag, scs_int max_iters,
                          ScsSolution *sol, ScsInfo *info, scs_int *rollbacks,
                          scs_int *mults_one, scs_int *found) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs_float *dr = SCS_NULL, *dr_init;
  ScsWork *w;
  scs_int status = SCS_FAILED, i, l;
  *rollbacks = *mults_one = *found = -1;
  memset(info, 0, sizeof(ScsInfo));
  g_acc_data(d, k, opt_sol);
  scs_set_default_settings(stgs);
  stgs->scale = G_RB_SCALE;
  stgs->adaptive_scale = adaptive;
  stgs->adaptive_diag_scale = diag;
  stgs->max_iters = max_iters;
  stgs->verbose = 0;
  w = scs_init(d, k, stgs);
  if (w) {
    l = d->n + d->m + 1;
    dr = (scs_float *)scs_calloc(l, sizeof(scs_float));
    memcpy(dr, w->diag_r, l * sizeof(scs_float));
    dr_init = w->diag_r;
    w->diag_r = dr;
    status = g_acc_solve_capture(
        w, sol, info, "WARNING: scale update at iteration",
        "rolled back: K g = (c, -b) could not be solved accurately with the "
        "new R",
        "it indicates a bug in SCS", found);
    *rollbacks = w->scale_rollbacks;
    *mults_one = 1;
    for (i = 0; w->scale_mults && i < d->m; ++i) {
      if (w->scale_mults[i] != 1.) {
        *mults_one = 0;
      }
    }
    w->diag_r = dr_init; /* the array scs_finish frees */
    scs_finish(w);
  }
  scs_free(dr);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  return status;
}

/* How far solve b's iterates are from solve a's: the largest over x, y and s
 * of ||a - b||_inf / max(||a||_inf, ||b||_inf) (0 where both are 0), in
 * double; NaN if an entry of either is NaN or infinite, so that the
 * comparison with G_RB_TOL fails. */
static double g_rb_dev(scs_int n, scs_int m, const ScsSolution *a,
                       const ScsSolution *b) {
  const scs_float *va[3], *vb[3];
  scs_int len[3], i, j;
  double dev = 0., num, den;
  va[0] = a->x, vb[0] = b->x, len[0] = n;
  va[1] = a->y, vb[1] = b->y, len[1] = m;
  va[2] = a->s, vb[2] = b->s, len[2] = m;
  for (j = 0; j < 3; ++j) {
    num = den = 0.;
    for (i = 0; i < len[j]; ++i) {
      if (!isfinite((double)va[j][i]) || !isfinite((double)vb[j][i])) {
        return NAN;
      }
      num = MAX(num, fabs((double)va[j][i] - (double)vb[j][i]));
      den = MAX(den, MAX(fabs((double)va[j][i]), fabs((double)vb[j][i])));
    }
    if (num > 0.) {
      dev = MAX(dev, num / den);
    }
  }
  return dev;
}

static const char *test_g_rollback(void) {
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *ref = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0}, info_ref = {0};
  scs_int run, status, status_ref, rollbacks, mults_one, found, rb_ref,
      mults_ref, found_ref, nnz, skip;
  double dev;
  const char *fail = SCS_NULL, *method = scs_get_lin_sys_method();

  /* the instance's dimensions, for the comparisons and the skip rule */
  g_acc_data(d, k, opt_sol);
  nnz = d->n + d->m + d->A->p[d->n] + d->P->p[d->n]; /* K's stored triangle */
  skip = (strstr(method, "cuDSS") && nnz >= CUDSS_MIN_PRODUCT_NNZ);
#if MKL_NO_OWN_PRODUCTS <= 0
  skip = skip || strstr(method, "mkl");
#endif
  if (skip) {
    scs_printf("test_g_rollback: %s measures with its own K, which follows "
               "the update: skipped\n",
               method);
  }
  for (run = 0; run < 2 && !skip && !fail; ++run) {
    /* run 0 (1. above): adaptive_diag_scale off; run 1 (2.): on */
    scs_int max_iters = run == 0 ? G_RB_ITERS : G_RB_ITERS_DIAG;
    status = g_rb_solve(1, run, max_iters, sol, &info, &rollbacks, &mults_one,
                        &found);
    status_ref = g_rb_solve(0, run, max_iters, ref, &info_ref, &rb_ref,
                            &mults_ref, &found_ref);
    dev = (status_ref == SCS_FAILED || status == SCS_FAILED)
              ? NAN
              : g_rb_dev(d->n, d->m, sol, ref);
    scs_printf("test_g_rollback: %s, adaptive_diag_scale %li: status %li "
               "after %li iterations, %li rollbacks, %li scale updates, "
               "scale %.2e; with adaptive_scale off: status %li after %li "
               "iterations; iterates apart by %.2e (%s)\n",
               method, (long)run, (long)status, (long)info.iter,
               (long)rollbacks, (long)info.scale_updates, (double)info.scale,
               (long)status_ref, (long)info_ref.iter, dev,
               run == 0 ? "not compared in this run" : "at most G_RB_TOL");
    if (info.iter != max_iters) {
      fail = "test_g_rollback: the solve did not run all its iterations (a "
             "failed restore stops it with SCS_FAILED)";
    } else if (rollbacks < 1) {
      fail = "test_g_rollback: no rollback, the test is vacuous";
    } else if (rollbacks * RESCALING_MIN_ITERS >= info.iter) {
      fail = "test_g_rollback: more than one rollback per "
             "RESCALING_MIN_ITERS iterations";
    } else if (info.scale_updates != 0 ||
               info.scale != (scs_float)G_RB_SCALE) {
      fail = "test_g_rollback: a scale update was kept";
    } else if (!mults_one) {
      fail = "test_g_rollback: the per-row multipliers were not put back";
    } else if (!found) {
      fail = "test_g_rollback: no WARNING line for the rollback in the "
             "output, or a message asking for a bug report";
    } else if (info_ref.iter != max_iters || rb_ref != 0) {
      fail = "test_g_rollback: the solve with adaptive_scale off failed";
    } else if (run == 1 && info.status_val != info_ref.status_val) {
      fail = "test_g_rollback: not the status of the solve with "
             "adaptive_scale off";
    } else if (run == 1 && !(dev <= G_RB_TOL)) {
      fail = "test_g_rollback: the iterates are more than G_RB_TOL from "
             "those of the solve with adaptive_scale off";
    }
  }
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(opt_sol);
  SCS(free_sol)(sol);
  SCS(free_sol)(ref);
  return fail;
}

/*
 * test_factor_rollback: a scale update whose factorization fails at the new
 * R is rolled back like one whose solve for g fails (scale_update_rollback in
 * src/scs.c). Any nonzero status of scs_update_lin_sys_diag_r is a failure
 * (linsys.h), so the dense solver's positive dpotrf info counts too. Before
 * this change (and in SCS 3.3.2) update_scale tested < 0: it stopped with
 * SCS_FAILED on QDLDL's -1 and ignored the dense solver's 1. On this problem
 * (measured on the Mac, before this change) the dense solver then kept two
 * updates whose factorizations had failed and iterated with the unfactorized
 * matrix, while QDLDL stopped.
 *
 * How the factorization is made to fail at a scale update and not at setup,
 * through the problem data alone. K = [R_x + P, A'; A, -R_y] has n positive
 * eigenvalues exactly when R_x + P + A' R_y^-1 A is positive definite, which
 * the direct solvers check at every factorization, at setup and at every
 * update (at setup: test_non_psd_p): QDLDL through its pivots, the dense
 * solver through its Cholesky factorization, MKL Pardiso and cuDSS through
 * the numbers of positive and negative pivots they report, Apple Accelerate
 * through the signs of its pivots (the diagonal of D). With a P that is
 * not positive semidefinite, that depends on R_y. Here P = diag(G_FR_P11, 1)
 * over the box -1e5 <= x <= 1e5 (A = G_FR_A [I; -I], b = 1, c = 1,
 * normalize off), so that matrix is diagonal, its first entry R_x + G_FR_P11
 * plus G_FR_A^2 R_y^-1 of the two rows of x_1. At setup every R_y^-1 is the
 * scale setting G_FR_SCALE: the entry is 2e-10 * 1e11 - 2 = 18 > 0, and K
 * factorizes. Any update from a scale above MAX_SCALE_VALUE (1e6) is clamped
 * to at most that, and the per-row multipliers to at most DIAG_SCALE_MULT_MAX
 * (1e3), so every R_y^-1 an update can install is at most 1e9 and the entry
 * at most 2e-10 * 1e9 - 2 = -1.8 < 0: every attempt fails, whatever the
 * iterates. The two R_y^-1 terms are 10 times |G_FR_P11| at setup and at
 * most a tenth of it at an attempt. The problem is not convex; the solve is a
 * fixture for the rollback, not a solve.
 *
 * Why the entries of A are small, rather than P11 large: setup needs
 * R_y < 1e-9 (above), and K's pivots must not be small next to its norm.
 * MKL Pardiso replaces a pivot below 10^-iparm[9] = 1e-13 times the infinity
 * norm of K by that size (with its sign), and cuDSS one below its pivot
 * epsilon (1e-13 in double precision): the factorization at setup would
 * then be that of another matrix. With A = [I; -I], P11 must be below -2e9,
 * so the norm is above 2e9 and Pardiso's threshold above 2e-4, far above
 * R_y. Here the norm is about 2 and every pivot of K at setup is at least
 * about 1e-11 in magnitude (R_y, or G_FR_A^2 / |G_FR_P11| = 5e-11, or sums
 * of them): 50 times Pardiso's threshold and 100 times cuDSS's epsilon.
 * cuDSS's epsilon in single precision (1e-5) is above every R_y that setup
 * can have here, so single-precision cuDSS cannot run this test as it is
 * (not built in CI).
 *
 * Checked, with each solver that refuses such a K at an update:
 * 1. G_FR_ITERS iterations with the default settings otherwise (verbose off):
 *    all iterations run (a failed restore stops the solve at once, with
 *    SCS_FAILED), at least one rollback, at most one per RESCALING_MIN_ITERS,
 *    no update kept, the scale still G_FR_SCALE, the per-row multipliers back
 *    at 1, a status other than SCS_FAILED, the WARNING line saying that the
 *    linear solver could not factorize K with the new R (printed whatever
 *    verbose is; the status it gives is the solver's: -1 for QDLDL, Apple
 *    Accelerate, MKL Pardiso and cuDSS, dpotrf's info, 1, for the dense
 *    solver), and no message asking for a bug report.
 * 2. The same with SCS's copy of rho_x set to G_FR_BAD_RHO_X after scs_init
 *    (written into the workspace, as test_g_rollback writes w->diag_r: setup
 *    would refuse it), so that both diagonal entries are negative at every R
 *    an update can install and at the previous one: the restore fails too,
 *    and SCS stops with SCS_FAILED and the ERROR line saying that the linear
 *    solver could not factorize K with the previous R either, without a
 *    rollback and without the message asking for a bug report.
 * The messages are checked unless printing is compiled out (NO_PRINTING).
 * The iterates are not compared with those of the same solve with
 * adaptive_scale off; test_g_rollback compares the continuation. (Measured on
 * the Mac, QDLDL and dense in double precision, after 150, 300 and 1000
 * iterations: the same bits; one ulp in one entry of b moved the solve with
 * adaptive_scale off by at most 8e-8, relative.)
 *
 * Skipped with the CPU and GPU CG solvers, with the reason printed: they
 * factorize nothing (their update builds a preconditioner and resets the
 * deflation basis, and fails only if the GPU one cannot copy R or the
 * preconditioner to the GPU; at an update only the check of the solve for g
 * applies to them, test_g_rollback). With CG (measured on the Mac) the first
 * update is kept and the solve continues with that K.
 */
/* G_FR_A, G_FR_P11, G_FR_SCALE: see above. G_FR_ITERS: room for 2 attempts,
 * at iterations 100 and 200 (measured: 2 rollbacks, QDLDL, dense and Apple
 * Accelerate, double and single precision, Mac). G_FR_BAD_RHO_X: 500 times
 * the largest G_FR_A^2 R_y^-1 term of setup (20 over both rows), so that
 * R_x + P + A' R_y^-1 A is negative definite at every R. */
#define G_FR_A (1e-5)
#define G_FR_P11 (-2.)
#define G_FR_SCALE (1e11)
#define G_FR_ITERS (300)
#define G_FR_BAD_RHO_X (-1e4)

/* One solve of the problem above (rho_x replaced by bad_rho_x after scs_init
 * when it is nonzero), with the output captured as in g_acc_solve_capture
 * (want, want2 and unwanted as there). Returns the status; *rollbacks is
 * w->scale_rollbacks, *mults_one whether every per-row multiplier is 1 at the
 * end (1 if there are none), or -1 for every output if setup failed. */
static scs_int g_fr_solve(scs_float bad_rho_x, ScsSolution *sol, ScsInfo *info,
                          const char *want, const char *want2,
                          const char *unwanted, scs_int *rollbacks,
                          scs_int *mults_one, scs_int *found) {
  scs_float Ax[] = {G_FR_A, -G_FR_A, G_FR_A, -G_FR_A};
  scs_int Ai[] = {0, 2, 1, 3};
  scs_int Ap[] = {0, 2, 4};
  scs_float Px[] = {G_FR_P11, 1.};
  scs_int Pi[] = {0, 1};
  scs_int Pp[] = {0, 1, 2};
  scs_float b[] = {1., 1., 1., 1.};
  scs_float c[] = {1., 1.};
  ScsMatrix A, P;
  ScsData d;
  ScsCone k;
  ScsSettings stgs;
  ScsWork *w;
  scs_int status = SCS_FAILED, i;
  *rollbacks = *mults_one = *found = -1;
  memset(info, 0, sizeof(ScsInfo));
  A.x = Ax, A.i = Ai, A.p = Ap, A.m = 4, A.n = 2;
  P.x = Px, P.i = Pi, P.p = Pp, P.m = 2, P.n = 2;
  d.m = 4, d.n = 2, d.A = &A, d.P = &P, d.b = b, d.c = c;
  memset(&k, 0, sizeof(ScsCone));
  k.l = 4;
  scs_set_default_settings(&stgs);
  stgs.normalize = 0;
  stgs.scale = G_FR_SCALE;
  stgs.max_iters = G_FR_ITERS;
  stgs.verbose = 0;
  w = scs_init(&d, &k, &stgs);
  if (w) {
    if (bad_rho_x != 0.) {
      w->stgs->rho_x = bad_rho_x;
    }
    status = g_acc_solve_capture(w, sol, info, want, want2, unwanted, found);
    *rollbacks = w->scale_rollbacks;
    *mults_one = 1;
    for (i = 0; w->scale_mults && i < d.m; ++i) {
      if (w->scale_mults[i] != 1.) {
        *mults_one = 0;
      }
    }
    scs_finish(w);
  }
  return status;
}

static const char *test_factor_rollback(void) {
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int status, rollbacks, mults_one, found;
  const char *fail = SCS_NULL, *method = scs_get_lin_sys_method();

  if (strstr(method, "indirect")) {
    scs_printf("test_factor_rollback: %s does not refuse a K with fewer than "
               "n positive eigenvalues at a scale update (the CG solvers "
               "factorize nothing): skipped\n",
               method);
    SCS(free_sol)(sol);
    return mu_skipped;
  }
  /* 1. every attempt fails to factorize and is rolled back */
  status = g_fr_solve(0., sol, &info, "WARNING: scale update at iteration",
                      "rolled back: the linear solver could not factorize K "
                      "with the new R",
                      "it indicates a bug in SCS", &rollbacks, &mults_one,
                      &found);
  scs_printf("test_factor_rollback: %s: status %li after %li iterations, %li "
             "rollbacks, %li scale updates, scale %.2e\n",
             method, (long)status, (long)info.iter, (long)rollbacks,
             (long)info.scale_updates, (double)info.scale);
  if (status == SCS_FAILED || info.status_val == SCS_FAILED) {
    fail = "test_factor_rollback: SCS_FAILED: the failed factorization was "
           "not rolled back";
  } else if (info.iter != G_FR_ITERS) {
    fail = "test_factor_rollback: the solve did not run all its iterations";
  } else if (rollbacks < 1) {
    fail = "test_factor_rollback: no rollback, the test is vacuous";
  } else if (rollbacks * RESCALING_MIN_ITERS >= info.iter) {
    fail = "test_factor_rollback: more than one rollback per "
           "RESCALING_MIN_ITERS iterations";
  } else if (info.scale_updates != 0 ||
             info.scale != (scs_float)G_FR_SCALE) {
    fail = "test_factor_rollback: a scale update was kept (a failed "
           "factorization not caught?)";
  } else if (!mults_one) {
    fail = "test_factor_rollback: the per-row multipliers were not put back";
  } else if (!found) {
    fail = "test_factor_rollback: no WARNING line saying that K could not be "
           "factorized with the new R, or a message asking for a bug report";
  }

  if (!fail) {
    /* 2. the restore fails too */
    status = g_fr_solve(G_FR_BAD_RHO_X, sol, &info,
                        "ERROR: scale update at iteration",
                        "could not factorize K with the previous R either",
                        "it indicates a bug in SCS", &rollbacks, &mults_one,
                        &found);
    scs_printf("test_factor_rollback: %s, rho_x %.0e after setup: status %li, "
               "%li rollbacks\n",
               method, G_FR_BAD_RHO_X, (long)status, (long)rollbacks);
    if (status != SCS_FAILED || info.status_val != SCS_FAILED) {
      fail = "test_factor_rollback: a failed restore did not stop with "
             "SCS_FAILED";
    } else if (rollbacks != 0) {
      fail = "test_factor_rollback: a failed restore was counted as a "
             "rollback";
    } else if (!found) {
      fail = "test_factor_rollback: a failed restore: no ERROR line saying "
             "that K could not be factorized with the previous R either, or a "
             "message asking for a bug report";
    }
  }
  SCS(free_sol)(sol);
  return fail;
}
