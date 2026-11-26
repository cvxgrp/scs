#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * The convergence check and the scale update read only iterates the Anderson
 * acceleration (AA) safeguard has vetted (see the AA safeguard in scs_solve).
 * When the safeguard rejects the AA step of an iteration with a check due,
 * the check runs at the next iteration instead, which never carries an AA
 * step; and no AA step is made at the last iteration.
 *
 * The problems are scaled-down copies of HUESTIS (Maros-Meszaros): minimize
 * ||x||^2 subject to two dense rows and x >= 0, where entry j of the rows is
 * the integral of t^2 (resp. t^4) over [j, j + 1], divided by n^3 (resp.
 * n^5), with HUESTIS's right-hand side. With n = 10,000 this is HUESTIS:
 * its entries run from 2e-21 to 1e-4, and the iterate first heads for an
 * infeasibility certificate while the scale climbs to its limit. Many AA
 * steps are proposed and rejected there, at checks too.
 *
 * Before the change, the check and update_scale read the iterate of the AA
 * step before the safeguard judged it; an update applied there reset AA, so
 * the safeguard did not test the step, and the step was kept. With QDLDL on
 * macOS arm64 (Accelerate for AA), n = 3,000, eps 1e-4 and 1e-6 and 20
 * initial scales each (the default moved by 0 to 19 ulps), 22 of the 40
 * solves did not solve within 22,600 iterations (54 to 131 scale updates),
 * and 14 of the 18 that did took 6,500 to 22,550 iterations with 36 to 90
 * updates. With the change all 40 solved in 3,725 to 5,250 iterations with
 * at most 12 updates. (HUESTIS itself: 71 of 82 runs failed with QDLDL on an
 * ARM Linux machine with NVPL; see scs_solve.)
 *
 * The tests read what happened at each iteration from the CSV log
 * (log_csv_filename): the scale and the counts of accepted and rejected AA
 * steps. Logging does not change the iterates (test_csv_log: the CG
 * tolerance, like the refinement's, is saved at the checks), so these are
 * the trajectories of unlogged solves too; and the structure asserted below
 * holds on any trajectory.
 */

/* the problem of size n described above */
static void aa_vetted_data(ScsData *d, ScsCone *k, scs_int n) {
  scs_int j, m = n + 3, nnz = 3 * n;
  scs_float n3 = (scs_float)n * n * n, n5 = n3 * n * n;
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsMatrix *P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  A->m = m;
  A->n = n;
  A->x = (scs_float *)scs_calloc(nnz, sizeof(scs_float));
  A->i = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  A->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  P->m = n;
  P->n = n;
  P->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  P->i = (scs_int *)scs_calloc(n, sizeof(scs_int));
  P->p = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  d->m = m;
  d->n = n;
  d->A = A;
  d->P = P;
  d->b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  /* rows 0, 1: the two equality rows; row 2: the box cone's t (= 1); rows 3
   * to n + 2: -x in the box [0, inf) */
  for (j = 0; j < n; ++j) {
    scs_float a = (scs_float)j, e = (scs_float)(j + 1);
    A->p[j] = 3 * j;
    A->x[3 * j] = (e * e * e - a * a * a) / 3. / n3;
    A->i[3 * j] = 0;
    A->x[3 * j + 1] = (e * e * e * e * e - a * a * a * a * a) / 5. / n5;
    A->i[3 * j + 1] = 1;
    A->x[3 * j + 2] = -1.;
    A->i[3 * j + 2] = 3 + j;
    P->x[j] = 2.;
    P->i[j] = j;
    P->p[j] = j;
  }
  A->p[n] = nnz;
  P->p[n] = n;
  d->b[0] = 1835.2;
  d->b[1] = 909.8;
  d->b[2] = 1.;
  k->z = 2;
  k->bsize = n + 1;
  k->bl = (scs_float *)scs_calloc(n, sizeof(scs_float));
  k->bu = (scs_float *)scs_calloc(n, sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    k->bu[j] = INFINITY;
  }
}

/* the settings of one solve: the defaults, eps_abs = eps_rel = eps, the
 * initial scale moved up by ulps ulps, and the AA settings given */
static void aa_vetted_settings(ScsSettings *stgs, scs_float eps, scs_int ulps,
                               scs_int max_iters, scs_int interval,
                               scs_int lookback, scs_int type1,
                               const char *csv) {
  scs_int u;
  scs_set_default_settings(stgs);
  stgs->eps_abs = eps;
  stgs->eps_rel = eps;
  stgs->max_iters = max_iters;
  stgs->acceleration_interval = interval;
  stgs->acceleration_lookback = lookback;
  stgs->acceleration_type_1 = type1;
  stgs->verbose = 0;
  stgs->log_csv_filename = csv;
  for (u = 0; u < ulps; ++u) {
#ifdef SFLOAT
    stgs->scale = nextafterf(stgs->scale, INFINITY);
#else
    stgs->scale = nextafter(stgs->scale, INFINITY);
#endif
  }
}

#if NO_READ_WRITE == 0
/* The CSV log of one solve: per row (iteration), the scale and the counts
 * of accepted and rejected AA steps so far. Row j is iteration j: the loop
 * logs iterations 0 to the last but one, and the last (the one the solve
 * stopped at, or max_iters) is logged after the loop. */
typedef struct {
  scs_int rows;
  scs_float *scale;
  scs_int *acc, *rej;
} AaVettedLog;

static void aa_vetted_log_free(AaVettedLog *g) {
  scs_free(g->scale);
  scs_free(g->acc);
  scs_free(g->rej);
  memset(g, 0, sizeof(AaVettedLog));
}

/* the start of field col of a CSV line (SCS_NULL if the line is shorter) */
static const char *aa_vetted_field(const char *line, scs_int col) {
  for (; col > 0 && line; --col) {
    line = strchr(line, ',');
    if (line) {
      line++;
    }
  }
  return line;
}

/* reads the log of a solve of at most max_iters iterations; returns 0, or -1
 * if the file is missing or malformed */
static scs_int aa_vetted_log_read(const char *csv, scs_int max_iters,
                                  AaVettedLog *g) {
  enum { LINE = 16384 };
  static char line[LINE];
  static const char *names[4] = {"iter", "scale", "accepted_accel_steps",
                                 "rejected_accel_steps"};
  scs_int col[4] = {-1, -1, -1, -1}, j, q, cap = max_iters + 1;
  const char *f;
  FILE *fin = fopen(csv, "r");
  memset(g, 0, sizeof(AaVettedLog));
  if (!fin) {
    return -1;
  }
  g->scale = (scs_float *)scs_calloc(cap, sizeof(scs_float));
  g->acc = (scs_int *)scs_calloc(cap, sizeof(scs_int));
  g->rej = (scs_int *)scs_calloc(cap, sizeof(scs_int));
  if (!fgets(line, LINE, fin)) {
    fclose(fin);
    return -1;
  }
  /* the header: one name per field, each followed by a comma */
  for (j = 0, f = line; f && *f && *f != '\n'; ++j, f = aa_vetted_field(f, 1)) {
    for (q = 0; q < 4; ++q) {
      size_t len = strlen(names[q]);
      if (!strncmp(f, names[q], len) && f[len] == ',') {
        col[q] = j;
      }
    }
  }
  for (q = 0; q < 4; ++q) {
    if (col[q] < 0) {
      fclose(fin);
      return -1;
    }
  }
  while (fgets(line, LINE, fin)) {
    const char *fi = aa_vetted_field(line, col[0]);
    const char *fs = aa_vetted_field(line, col[1]);
    const char *fa = aa_vetted_field(line, col[2]);
    const char *fr = aa_vetted_field(line, col[3]);
    if (g->rows >= cap || !fi || !fs || !fa || !fr ||
        strtol(fi, SCS_NULL, 10) != g->rows) {
      fclose(fin);
      return -1;
    }
    g->scale[g->rows] = strtod(fs, SCS_NULL);
    g->acc[g->rows] = (scs_int)strtod(fa, SCS_NULL);
    g->rej[g->rows] = (scs_int)strtod(fr, SCS_NULL);
    g->rows++;
  }
  fclose(fin);
  return g->rows > 0 ? 0 : -1;
}

/* steps accepted (rejected) at iteration j */
static scs_int aa_vetted_acc(const AaVettedLog *g, scs_int j) {
  return j <= 0 ? 0 : g->acc[j] - g->acc[j - 1];
}
static scs_int aa_vetted_rej(const AaVettedLog *g, scs_int j) {
  return j <= 0 ? 0 : g->rej[j] - g->rej[j - 1];
}

/* what the tests count over many solves */
typedef struct {
  scs_int solves, changes;
  /* scale changes at an iteration that is not a multiple of
   * CONVERGED_INTERVAL: checks that waited for a rejected step */
  scs_int waited;
  /* scale changes at an iteration whose AA step was accepted */
  scs_int on_accepted;
  /* solves that stopped at a check that waited */
  scs_int stops_waited;
} AaVettedCount;

/* Checks the structure of one solve's log (returns SCS_NULL, or what does
 * not hold), and adds to the counts. It holds on any trajectory:
 * (a) the scale changes only at a check, which is at a multiple of
 *     CONVERGED_INTERVAL whose AA step was not rejected, or one iteration
 *     after a multiple whose step was;
 * (b) the iteration after a rejected step carries no step;
 * (c) a solve that converged stopped at such a check;
 * (d) no step is made at the last iteration, max_iters - 1. */
static const char *aa_vetted_structure(const AaVettedLog *g,
                                       const ScsInfo *info, scs_int max_iters,
                                       AaVettedCount *c) {
  scs_int j, last = g->rows - 1;
  c->solves++;
  if (last != info->iter) {
    return "the log does not end at the iteration the solve stopped at";
  }
  for (j = 1; j <= last; ++j) {
    if (aa_vetted_rej(g, j - 1) &&
        (aa_vetted_rej(g, j) || aa_vetted_acc(g, j))) {
      return "the iteration after a rejected AA step carried a step";
    }
    if (g->scale[j] == g->scale[j - 1]) {
      continue;
    }
    c->changes++;
    if (aa_vetted_rej(g, j)) {
      return "the scale changed at an iteration whose AA step was rejected";
    }
    if (j % CONVERGED_INTERVAL != 0) {
      if (!((j - 1) % CONVERGED_INTERVAL == 0 && aa_vetted_rej(g, j - 1))) {
        return "the scale changed at an iteration with no check";
      }
      c->waited++;
    }
    c->on_accepted += aa_vetted_acc(g, j);
  }
  if (info->status_val == SCS_SOLVED || info->status_val == SCS_INFEASIBLE ||
      info->status_val == SCS_UNBOUNDED) {
    if (aa_vetted_rej(g, last)) {
      return "the solve stopped at an iteration whose AA step was rejected";
    }
    if (last % CONVERGED_INTERVAL != 0) {
      if (!((last - 1) % CONVERGED_INTERVAL == 0 &&
            aa_vetted_rej(g, last - 1))) {
        return "the solve stopped at an iteration with no check";
      }
      c->stops_waited++;
    }
  }
  if (info->iter == max_iters && last >= 2 &&
      (aa_vetted_acc(g, max_iters - 1) || aa_vetted_rej(g, max_iters - 1))) {
    return "an AA step was made at the last iteration";
  }
  return SCS_NULL;
}

#define AA_VETTED_CSV "test_aa_vetted.csv"

/* one logged solve of the problem of size n; checks its structure and adds
 * to the counts */
static const char *aa_vetted_logged(scs_int n, scs_float eps, scs_int ulps,
                                    scs_int max_iters, scs_int interval,
                                    scs_int lookback, scs_int type1,
                                    AaVettedCount *c, ScsInfo *info) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  AaVettedLog g;
  const char *fail = SCS_NULL;
  memset(info, 0, sizeof(ScsInfo));
  aa_vetted_data(d, k, n);
  aa_vetted_settings(stgs, eps, ulps, max_iters, interval, lookback, type1,
                     AA_VETTED_CSV);
  scs(d, k, stgs, sol, info);
  if (aa_vetted_log_read(AA_VETTED_CSV, max_iters, &g) != 0) {
    fail = "could not read the CSV log";
  } else {
    fail = aa_vetted_structure(&g, info, max_iters, c);
  }
  if (fail) {
    scs_printf("test_aa_vetted: n %li, eps %.0e, interval %li, lookback %li, "
               "type %s, max_iters %li: %s\n",
               (long)n, eps, (long)interval, (long)lookback,
               type1 ? "I" : "II", (long)max_iters, fail);
  }
  aa_vetted_log_free(&g);
  remove(AA_VETTED_CSV);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  return fail;
}
#endif

/*
 * The structure (aa_vetted_structure) on every solve, and both sides of the
 * safeguard's verdict in these solves:
 * - checks that waited for a rejected step and changed the scale there
 *   (before the change a check never waited);
 * - scale changes at iterations whose AA step was accepted, and at
 *   acceleration_interval 1 too, where every check is an AA iteration.
 *   (Skipping the scale update at every iteration that carries an AA step
 *   instead, whatever the safeguard says, would never update the scale at
 *   such an iteration, and at interval 1 would update it rarely while AA's
 *   steps are being accepted.)
 * The solves: at interval 1, n = 300, 500 and 700, lookbacks 1 and 2 (where
 * AA's steps are accepted most often on these problems), both AA types and two
 * initial scales (24 solves); at intervals 2, 3 and 5 (the default), n = 300,
 * lookbacks 1, 3 and 10 (the default) and both AA types (18 solves); eps 1e-4
 * and at most 2,000 iterations each. Measured on macOS arm64 with QDLDL (and
 * CPU CG): at interval 1, 34 (32) checks that waited and changed the scale
 * and 22 (16) scale changes at an accepted step; at the other intervals, 15
 * (15) and 12 (12). The test asks for at least one of each, and at interval 1
 * for at least one change at an accepted step. Before the change no check
 * ever waited (and this test stops earlier, at an AA step made at the last
 * iteration of a solve that ran out of iterations).
 * Without LAPACK there is no AA, and with CONVERGED_INTERVAL 1 (VERBOSITY >
 * 0) every iteration is a check, so then only the structure is checked.
 */
static const char *test_aa_vetted_structure(void) {
#if NO_READ_WRITE > 0
  return mu_skipped;
#else
  static const scs_int n1[] = {300, 500, 700}, intervals[] = {2, 3, 5},
                       lookbacks[] = {1, 3, 10};
  AaVettedCount one = {0}, all = {0};
  ScsInfo info;
  scs_int a, b, t, u;
  const char *fail;
  for (a = 0; a < 3; ++a) {
    for (b = 1; b <= 2; ++b) {
      for (t = 0; t < 2; ++t) {
        for (u = 0; u < 2; ++u) {
          fail = aa_vetted_logged(n1[a], 1e-4, u, 2000, 1, b, t, &one, &info);
          if (fail) {
            return fail;
          }
        }
      }
    }
  }
  for (a = 0; a < 3; ++a) {
    for (b = 0; b < 3; ++b) {
      for (t = 0; t < 2; ++t) {
        fail = aa_vetted_logged(300, 1e-4, 0, 2000, intervals[a], lookbacks[b],
                                t, &all, &info);
        if (fail) {
          return fail;
        }
      }
    }
  }
  scs_printf("test_aa_vetted_structure: interval 1: %li solves, %li scale "
             "changes, %li after a wait, %li at an accepted step, %li stops "
             "after a wait; intervals 2, 3, 5: %li solves, %li changes, %li "
             "after a wait, %li at an accepted step, %li stops after a wait\n",
             (long)one.solves, (long)one.changes, (long)one.waited,
             (long)one.on_accepted, (long)one.stops_waited, (long)all.solves,
             (long)all.changes, (long)all.waited, (long)all.on_accepted,
             (long)all.stops_waited);
#if defined(USE_LAPACK) && CONVERGED_INTERVAL > 1
  mu_assert("test_aa_vetted_structure: no check waited for a rejected AA "
            "step and changed the scale",
            one.waited + all.waited > 0);
  mu_assert("test_aa_vetted_structure: the scale never changed at an "
            "iteration whose AA step was accepted",
            one.on_accepted + all.on_accepted > 0);
  mu_assert("test_aa_vetted_structure: at acceleration_interval 1 the scale "
            "never changed at an iteration whose AA step was accepted",
            one.on_accepted > 0);
#endif
  return SCS_NULL;
#endif
}

/*
 * No AA step at the last iteration, and the structure, when the solve runs
 * out of iterations with max_iters - 1 a multiple of the AA interval: n =
 * 1,000 at eps 1e-6 (which needs more iterations than any of these), at
 * interval 1 with lookbacks 1 (type II) and 2 (type I) and max_iters 101,
 * 151, ..., 401, and at max_iters 1,001 and 1,251 (interval 5) and 301
 * (interval 3) with the default lookback. Before the change, AA made a step
 * at the last iteration in 9 of these 17 solves with QDLDL and in 7 with CPU
 * CG (macOS arm64). A solve that converges first does not test the last
 * iteration (with CG, 5 of the 17 solved in 200 iterations), so the test asks
 * only that some solve ran out of iterations (all 17 with QDLDL, 12 with CG).
 */
static const char *test_aa_vetted_last_iter(void) {
#if NO_READ_WRITE > 0
  return mu_skipped;
#else
  static const scs_int cfg[][4] = {
      /* interval, lookback, type I, max_iters */
      {5, 10, 1, 1001}, {5, 10, 1, 1251}, {3, 10, 1, 301}};
  AaVettedCount c = {0};
  ScsInfo info;
  scs_int j, b, mi, ran_out = 0;
  const char *fail;
  for (j = 0; j < 3 + 2 * 7; ++j) {
    scs_int iv = 1, lb, t1, max_iters;
    if (j < 3) {
      iv = cfg[j][0];
      lb = cfg[j][1];
      t1 = cfg[j][2];
      max_iters = cfg[j][3];
    } else {
      b = (j - 3) / 7;
      mi = (j - 3) % 7;
      lb = b + 1;
      t1 = b;
      max_iters = 101 + 50 * mi;
    }
    fail = aa_vetted_logged(1000, 1e-6, 0, max_iters, iv, lb, t1, &c, &info);
    if (fail) {
      return fail;
    }
    ran_out += info.iter == max_iters;
  }
  scs_printf("test_aa_vetted_last_iter: %li of 17 solves ran out of "
             "iterations\n",
             (long)ran_out);
  mu_assert("test_aa_vetted_last_iter: no solve ran out of iterations",
            ran_out > 0);
  return SCS_NULL;
#endif
}

/*
 * Two solves on one workspace, the second after scs_update has changed b:
 * the structure holds in each (a check that waits does not carry over to
 * the next solve). Whether they solve is printed, not checked: it is not what
 * this test is about, and it depends on the BLAS and LAPACK. At
 * acceleration_interval 1 this problem can take more than these 5,000
 * iterations: with the dense solver the CI runners (MKL on Linux, Accelerate
 * on macOS arm64) did not solve it from the default start, while on an Apple
 * M-series Mac every backend solved both solves from 41 starts, with
 * Accelerate and with OpenBLAS (the first in at most 3,551 iterations). In
 * single precision it does not solve at eps 1e-4.
 */
static const char *test_aa_vetted_repeat(void) {
#if NO_READ_WRITE > 0
  return mu_skipped;
#else
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs_float *b_new;
  AaVettedCount c = {0};
  AaVettedLog g;
  ScsInfo info = {0};
  ScsWork *w;
  scs_int r, status[2] = {SCS_FAILED, SCS_FAILED}, iters[2] = {0, 0};
  const char *fail = SCS_NULL;
  aa_vetted_data(d, k, 1000);
  aa_vetted_settings(stgs, 1e-4, 0, 5000, 1, 10, 1, AA_VETTED_CSV);
  b_new = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  memcpy(b_new, d->b, d->m * sizeof(scs_float));
  b_new[0] *= 1.1;
  b_new[1] *= 0.9;
  w = scs_init(d, k, stgs);
  for (r = 0; w && r < 2 && !fail; ++r) {
    if (r == 1 && scs_update(w, b_new, SCS_NULL) != 0) {
      fail = "test_aa_vetted_repeat: scs_update failed";
      break;
    }
    status[r] = scs_solve(w, sol, &info, r);
    iters[r] = info.iter;
    if (aa_vetted_log_read(AA_VETTED_CSV, stgs->max_iters, &g) != 0) {
      fail = "test_aa_vetted_repeat: could not read the CSV log";
    } else {
      fail = aa_vetted_structure(&g, &info, stgs->max_iters, &c);
    }
    aa_vetted_log_free(&g);
    remove(AA_VETTED_CSV);
  }
  if (w) {
    scs_finish(w);
  } else {
    fail = "test_aa_vetted_repeat: scs_init failed";
  }
  scs_free(b_new);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  if (fail) {
    return fail;
  }
  scs_printf("test_aa_vetted_repeat: exit flags %li and %li, %li and %li "
             "iterations\n",
             (long)status[0], (long)status[1], (long)iters[0], (long)iters[1]);
  (void)status; /* NO_PRINTING */
  (void)iters;
  return SCS_NULL;
#endif
}

/*
 * The failure itself, with QDLDL: the problem of size 3,000 at eps 1e-4 and
 * 1e-6 from two initial scales (the default and one ulp above) must solve
 * within 15,000 iterations and with at most 30 scale updates, about three
 * times what the change needed on the 40 starts measured (at most 5,250
 * iterations and 12 updates; see the top of this file). Before the change
 * these four solves took 17,475 to 17,775 iterations with 66 to 77 updates
 * (three) or did not solve within 22,600 (one). Other linear solvers are not
 * asserted: CPU CG solved 39 of the 40 before the change (and all 40 with it,
 * in 2,800 to 4,626 iterations), and the others were not measured. The tiny
 * entries need double precision.
 */
static const char *test_aa_vetted_huestis(void) {
#ifdef SFLOAT
  return mu_skipped;
#else
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info;
  scs_int u, e, status, worst = 0, updates = 0;
  static const scs_float eps[] = {1e-4, 1e-6};
  if (!strstr(scs_get_lin_sys_method(), "qdldl")) {
    return mu_skipped;
  }
  for (e = 0; e < 2; ++e) {
    for (u = 0; u < 2; ++u) {
      k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
      d = (ScsData *)scs_calloc(1, sizeof(ScsData));
      stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
      sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
      memset(&info, 0, sizeof(ScsInfo));
      aa_vetted_data(d, k, 3000);
      aa_vetted_settings(stgs, eps[e], u, 15000, 5, 10, 1, SCS_NULL);
      status = scs(d, k, stgs, sol, &info);
      SCS(free_data)(d);
      SCS(free_cone)(k);
      SCS(free_sol)(sol);
      scs_free(stgs);
      scs_printf("test_aa_vetted_huestis: eps %.0e, %li ulps: %s, %li "
                 "iterations, %li scale updates\n",
                 eps[e], (long)u, info.status, (long)info.iter,
                 (long)info.scale_updates);
      mu_assert("test_aa_vetted_huestis: not solved", status == SCS_SOLVED);
      mu_assert("test_aa_vetted_huestis: more than 30 scale updates",
                info.scale_updates <= 30);
      worst = MAX(worst, info.iter);
      updates = MAX(updates, info.scale_updates);
    }
  }
  scs_printf("test_aa_vetted_huestis: at most %li iterations and %li scale "
             "updates\n",
             (long)worst, (long)updates);
  return SCS_NULL;
#endif
}
