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
 * solves for g to about machine precision: a residual ||[c; -b] - K g||_inf
 * of at most G_TOL_EPS machine epsilons of the size of the terms it is
 * computed from, || |[c; -b]| + |K| |g| ||_inf, and stops with SCS_FAILED,
 * asking for a bug report, if that solve fails (glbopts.h, "Accuracy of the
 * solve for g"). This test checks both, for every linear solver the test
 * binary is built with, on a badly scaled QP that is not normalized (so that
 * K keeps the bad scaling), with that measure computed here:
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

/* SCS refines g until its own residual is at most G_TOL_EPS machine epsilons
 * of the terms; the residual recomputed here rounds differently, by up to a
 * few epsilons of the terms, hence the factor of 4. A plain solve is far
 * above it (printed). */
#define G_ACC_LIMIT (4 * G_TOL_EPS * SCS_MACH_EPS)
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
 * computed here, in double whatever scs_float is. P is diagonal in this
 * instance. */
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
    for (p = P->p[j]; p < P->p[j + 1]; ++p) { /* diagonal */
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
 * want and, unless it is SCS_NULL, want2 (1 when printing is compiled out,
 * or if the capture cannot be set up: then the output goes to stdout as
 * usual and only the status is checked) */
static scs_int g_acc_solve_capture(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                                   const char *want, const char *want2,
                                   scs_int *found) {
  scs_int status;
#if NO_PRINTING > 0
  *found = 1;
  return scs_solve(w, sol, info, 0);
#else
  FILE *tmp = tmpfile();
  int saved = -1;
  char buf[4096];
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
           (!want2 || strstr(buf, want2) != SCS_NULL);
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
                                 SCS_NULL, &found);
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
          "open an issue at https://github.com/cvxgrp/scs/issues", &found);
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
            "open an issue at https://github.com/cvxgrp/scs/issues", &found);
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
