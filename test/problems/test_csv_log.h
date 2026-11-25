#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"
#include <stdio.h>
#include <string.h>

/*
 * Logging to CSV (log_csv_filename) recomputes SCS's residuals every
 * iteration, to log them, and must not change the iterates: the algorithm
 * reads only the residuals computed at its convergence checks and at the end
 * of a solve. Before this change SCS set each iteration's CG tolerance
 * (project_lin_sys in src/scs.c) from the residuals as they stood, so with
 * logging from the previous iteration's instead of the last check's. On the
 * instance below, with CG on Apple silicon, both logged solves ended at other
 * x, y and s than the unlogged ones (in single precision after 250 and 25
 * iterations instead of 325 and 150). Now the check saves the size the
 * tolerance is set from (set_cg_tol_res), as it saves the tol of the
 * per-iteration refinement (set_refine_tol).
 *
 * test_csv_log solves a random SOCP (zero, linear and second-order cone
 * rows), with Anderson acceleration and adaptive scale on (the defaults),
 * with and without logging, twice on each workspace: a cold solve, then
 * another cold solve after scs_update changes b, whose first iteration runs
 * before that solve's first check (its CG tolerance comes from the residuals
 * the first solve ended with). It requires the same iterations, status, scale
 * updates and solution, bit for bit.
 *
 * Bit for bit needs solves that give the same bits whenever they compute with
 * the same values, in whichever workspace. SCS's own linear solvers (CG and
 * QDLDL) do, given a BLAS that does not round the same values differently at
 * another address, so the test runs with those two and checks the BLAS first
 * (csv_log_blas_by_placement). MKL's BLAS can round by alignment (it did in
 * CI, see test_solve_refine_csv), and the two workspaces' arrays sit at
 * different addresses (logging copies the file name, one more allocation), so
 * with it two solves can differ in the last bits with or without logging, and
 * the test is skipped. The other linear solvers' libraries do not promise the
 * same bits from one call to the next (GPU kernels, MKL Pardiso, threaded
 * factorizations: Apple's Accelerate solver gave different results from one
 * run of this suite to the next), so the test is skipped with them too; the
 * code it checks, in src/scs.c, is the same for every linear solver.
 */

/* the instance: n variables, m = 4 n rows (10% zero cone, 30% linear, the
 * rest in second-order cones of sizes drawn below), 10 entries per column of
 * A, primal and dual feasible (gen_random_prob_data) */
#define CSV_LOG_N (100)
#define CSV_LOG_SEED (2)
#define CSV_LOG_FILE "test_csv_log.csv"

static void csv_log_data(ScsData **dp, ScsCone **kp) {
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsSolution *opt = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs_int n = CSV_LOG_N, m = 4 * CSV_LOG_N, col_nnz = 10, rows, size;
  scs_int max_q = 67; /* ceil(m / log(m)) */
  d->n = n;
  d->m = m;
  k->z = m / 10;
  k->l = (3 * m) / 10;
  rows = m - k->z - k->l;
  k->q = (scs_int *)scs_calloc(rows, sizeof(scs_int));
  ran_start(CSV_LOG_SEED + 1000);
  while (rows > max_q) {
    size = (scs_int)(ran_arr_next() % max_q) + 1;
    k->q[k->qsize++] = size;
    rows -= size;
  }
  if (rows > 0) {
    k->q[k->qsize++] = rows;
  }
  gen_random_prob_data(n * col_nnz, col_nnz, d, k, opt, CSV_LOG_SEED);
  SCS(free_sol)(opt);
  *dp = d;
  *kp = k;
}

/* Whether SCS's vector operations can round the same values differently at
 * another address, which an exact comparison of two workspaces cannot allow:
 * malloc promises their arrays 16-byte alignment, no more, and they sit at
 * different addresses. With USE_LAPACK in double precision SCS computes dot
 * products, 2-norms and axpy with the BLAS, and some BLAS split a sum by
 * address: MKL's by alignment, Apple's Accelerate in its dot products of 900
 * entries and more (on Apple silicon; not in those of 800 or fewer). So the
 * probe computes, from the same values at 16 offsets 16 bytes apart, SCS's
 * dot product and 2-norm at every length up to len (the length of SCS's
 * iterate) and its axpy at len, and compares each with the first offset's bit
 * for bit. Returns 1 if any differs (or the probe cannot run). It calls only
 * SCS's own functions (the CMake test targets do not link the BLAS), so it
 * does not see Anderson acceleration's matrix products, which go to the BLAS
 * in every precision, or eigCG deflation's (CG in double precision with
 * USE_LAPACK: gemm and LAPACK's syev); a BLAS that rounds those by address
 * but not its dot products would go unnoticed. Without USE_LAPACK it is all
 * SCS's own loops, and there is no acceleration or deflation. */
static scs_int csv_log_blas_by_placement(scs_int len) {
#ifndef USE_LAPACK
  (void)len;
  return 0;
#else
  const scs_int e = 16 / sizeof(scs_float), shifts = 16, pad = shifts * e;
  scs_int nr = 3 * len, i, sh, differ = 0;
  scs_float *a = (scs_float *)scs_calloc(len + pad, sizeof(scs_float));
  scs_float *b = (scs_float *)scs_calloc(len + pad, sizeof(scs_float));
  scs_float *va = (scs_float *)scs_calloc(len, sizeof(scs_float));
  scs_float *vb = (scs_float *)scs_calloc(len, sizeof(scs_float));
  scs_float *r = (scs_float *)scs_calloc(nr, sizeof(scs_float));
  scs_float *r0 = (scs_float *)scs_calloc(nr, sizeof(scs_float));
  if (!a || !b || !va || !vb || !r || !r0) {
    differ = 1;
    goto out;
  }
  ran_start(CSV_LOG_SEED);
  for (i = 0; i < len; ++i) {
    va[i] = rand_scs_float();
    vb[i] = rand_scs_float();
  }
  for (sh = 0; sh < shifts && !differ; ++sh) {
    scs_float *x = a + sh * e, *y = b + sh * e;
    memcpy(x, va, len * sizeof(scs_float));
    memcpy(y, vb, len * sizeof(scs_float));
    for (i = 1; i <= len; ++i) {
      r[2 * i - 2] = SCS(dot)(x, y, i);
      r[2 * i - 1] = SCS(norm_2)(y, i);
    }
    SCS(add_scaled_array)(y, x, len, 0.3); /* y += 0.3 x */
    memcpy(r + 2 * len, y, len * sizeof(scs_float));
    if (sh == 0) {
      memcpy(r0, r, nr * sizeof(scs_float));
    } else {
      differ = memcmp(r, r0, nr * sizeof(scs_float)) != 0;
    }
  }
out:
  scs_free(a);
  scs_free(b);
  scs_free(va);
  scs_free(vb);
  scs_free(r);
  scs_free(r0);
  return differ;
#endif
}

/* solves the instance twice on one workspace (see above), logging to csv
 * unless it is SCS_NULL: status, info and solution of each solve in
 * status[i], info[i], sol[i]. Returns 0, or -1 if the setup failed. */
static scs_int csv_log_run(const char *csv, scs_int *status, ScsInfo *info,
                           ScsSolution *sol) {
  ScsData *d;
  ScsCone *k;
  ScsSettings stgs;
  ScsWork *w;
  scs_float *b2;
  scs_int i;
  csv_log_data(&d, &k);
  scs_set_default_settings(&stgs);
  stgs.verbose = 0;
  /* the defaults, but the test is about them: update_scale reads the
   * residuals, and Anderson acceleration mixes in every earlier iterate */
  stgs.adaptive_scale = 1;
  stgs.acceleration_lookback = ACCELERATION_LOOKBACK;
  stgs.log_csv_filename = csv;
  w = scs_init(d, k, &stgs);
  b2 = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  if (!w || !b2) {
    scs_finish(w);
    scs_free(b2);
    SCS(free_data)(d);
    SCS(free_cone)(k);
    return -1;
  }
  status[0] = scs_solve(w, &sol[0], &info[0], 0);
  for (i = 0; i < d->m; ++i) {
    b2[i] = 1.01 * d->b[i];
  }
  scs_update(w, b2, SCS_NULL);
  status[1] = scs_solve(w, &sol[1], &info[1], 0);
  scs_finish(w);
  scs_free(b2);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  return 0;
}

static const char *test_csv_log(void) {
#if NO_READ_WRITE > 0
  return mu_skipped; /* logging is compiled out */
#else
  const char *method = scs_get_lin_sys_method(), *fail = SCS_NULL;
  scs_int n = CSV_LOG_N, m = 4 * CSV_LOG_N, status0[2], status1[2], j;
  ScsInfo info0[2], info1[2];
  ScsSolution sol0[2], sol1[2];
  FILE *f;
  long size = 0;
  if (!strstr(method, "sparse-indirect-scs") && !strstr(method, "qdldl")) {
    scs_printf("test_csv_log: skipped with %s, whose library does not promise "
               "the same bits from one solve to the next\n",
               method);
    return mu_skipped;
  }
  if (csv_log_blas_by_placement(n + m + 1)) {
    scs_printf("test_csv_log: skipped, this BLAS rounds the same values "
               "differently at another place in memory, so two workspaces "
               "can differ in the last bits with or without logging\n");
    return mu_skipped;
  }
  memset(sol0, 0, sizeof(sol0));
  memset(sol1, 0, sizeof(sol1));
  remove(CSV_LOG_FILE);
  if (csv_log_run(SCS_NULL, status0, info0, sol0) != 0 ||
      csv_log_run(CSV_LOG_FILE, status1, info1, sol1) != 0) {
    fail = "test_csv_log: setup failed";
    goto out;
  }
  f = fopen(CSV_LOG_FILE, "r");
  if (f) {
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fclose(f);
  }
  remove(CSV_LOG_FILE);
  for (j = 0; j < 2; ++j) {
    scs_printf("test_csv_log: %s, solve %li: without logging %li iterations "
               "(%s), with logging %li iterations (%s)\n",
               method, (long)j + 1, (long)info0[j].iter, info0[j].status,
               (long)info1[j].iter, info1[j].status);
  }
  if (size <= 0) {
    fail = "test_csv_log: nothing was logged";
  } else if (status0[0] != SCS_SOLVED || status0[1] != SCS_SOLVED) {
    fail = "test_csv_log: not solved";
  } else if (info0[0].iter <= CONVERGED_INTERVAL) {
    /* nothing would run between two checks */
    fail = "test_csv_log: solved at the first check";
  }
  for (j = 0; j < 2 && !fail; ++j) {
    if (status1[j] != status0[j] || info1[j].iter != info0[j].iter ||
        info1[j].scale_updates != info0[j].scale_updates ||
        memcmp(sol1[j].x, sol0[j].x, n * sizeof(scs_float)) != 0 ||
        memcmp(sol1[j].y, sol0[j].y, m * sizeof(scs_float)) != 0 ||
        memcmp(sol1[j].s, sol0[j].s, m * sizeof(scs_float)) != 0) {
      fail = j == 0 ? "test_csv_log: logging changed the first solve"
                    : "test_csv_log: logging changed the second solve";
    }
  }
out:
  for (j = 0; j < 2; ++j) {
    scs_free(sol0[j].x);
    scs_free(sol0[j].y);
    scs_free(sol0[j].s);
    scs_free(sol1[j].x);
    scs_free(sol1[j].y);
    scs_free(sol1[j].s);
  }
  return fail;
#endif
}
