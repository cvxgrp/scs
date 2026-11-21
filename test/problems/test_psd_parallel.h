#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* A build with USE_PARALLEL_PSD and OpenMP projects the PSD blocks in parallel
 * (cones.c). Each block goes through the same LAPACK calls with the same
 * work-array sizes as in the serial projection, so with single-threaded BLAS
 * the result is the same bit for bit, whatever the number of threads.
 * test_psd_parallel projects one vector with 1 and with 4 threads, and
 * test_psd_parallel_solve solves one problem both ways.
 *
 * Every other build projects the blocks one after another. With OpenMP and
 * without USE_PARALLEL_PSD, test_psd_parallel checks that the projection stays
 * serial under 4 OpenMP threads.
 *
 * The comparisons are bit for bit when the BLAS is an OpenBLAS that the tests
 * can set to one thread. They are to rounding error for the other libraries,
 * which may round one block differently from call to call (MKL: the result
 * can depend on the alignment of the workspace) or cannot be set to one
 * thread from here (Accelerate). */

#if defined(USE_PARALLEL_PSD) && defined(_OPENMP)
/* this build projects the PSD blocks in parallel */
#define PSD_PAR_BUILD (1)
#endif

#ifdef PSD_PAR_BUILD
#if defined(__GNUC__) && defined(__ELF__)
/* weak references, NULL when the library is not linked (ELF only) */
#define PSD_PAR_WEAK_BLAS (1)
#ifdef __cplusplus
extern "C" {
#endif
void openblas_set_num_threads(int) __attribute__((weak));
int openblas_get_num_threads(void) __attribute__((weak));
int openblas_get_parallel(void) __attribute__((weak));
void MKL_Set_Num_Threads(int) __attribute__((weak));
int MKL_Get_Max_Threads(void) __attribute__((weak));
int MKL_Set_Num_Threads_Local(int) __attribute__((weak));
#ifdef __cplusplus
}
#endif
#else
#define PSD_PAR_WEAK_BLAS (0) /* e.g. Accelerate: VECLIB_MAXIMUM_THREADS=1 */
#endif

/* the thread settings before a test, and what the test can rely on */
typedef struct {
  int omp_nt;          /* OpenMP threads */
  int ob_nt;           /* OpenBLAS threads, 0: no OpenBLAS found */
  int mkl_nt;          /* MKL threads, 0: no MKL found */
  scs_int exact;       /* projections are the same bit for bit */
  scs_int exact_solve; /* and so is everything else in a solve */
} PsdParState;

/* Save the thread settings in st. Then BLAS on one thread where the library
 * can be told: a threaded BLAS rounds differently in the serial projection,
 * and a pthreads OpenBLAS with more than one thread keeps the projection
 * serial. An OpenMP OpenBLAS follows the OpenMP thread count instead. */
static void psd_par_enter(PsdParState *st) {
  memset(st, 0, sizeof(PsdParState));
  st->omp_nt = omp_get_max_threads();
#if PSD_PAR_WEAK_BLAS
  if (openblas_set_num_threads && openblas_get_num_threads) {
    /* an OpenMP OpenBLAS (2) follows the OpenMP thread count, which the tests
     * set, in the BLAS calls outside of the PSD threads. It is not given a
     * count of its own: from 0.3.34 on it would keep that count and ignore
     * the OpenMP one, and test_psd_serial_hold could not see the hold. */
    scs_int omp_ob = openblas_get_parallel && openblas_get_parallel() == 2;
    if (!omp_ob) {
      st->ob_nt = openblas_get_num_threads();
      openblas_set_num_threads(1);
    }
    st->exact = !MKL_Set_Num_Threads;
    st->exact_solve = st->exact && openblas_get_parallel && !omp_ob;
  }
  if (MKL_Set_Num_Threads && MKL_Get_Max_Threads) {
    st->mkl_nt = MKL_Get_Max_Threads();
    MKL_Set_Num_Threads(1);
  }
#endif
}

/* back to the settings saved in st */
static void psd_par_leave(PsdParState *st) {
  omp_set_num_threads(st->omp_nt);
#if PSD_PAR_WEAK_BLAS
  if (st->ob_nt > 0)
    openblas_set_num_threads(st->ob_nt);
  if (st->mkl_nt > 0)
    MKL_Set_Num_Threads(st->mkl_nt);
#endif
}

/* Whether the blocks stay serial even when 4 OpenMP threads are asked for, so
 * that there is nothing to test (says why). Leaves the thread count at 4. */
static scs_int psd_par_unavailable(void) {
  omp_set_num_threads(4);
#ifdef PSD_PAR_MIN_COST
  scs_printf("skipping: PSD_PAR_MIN_COST is set in this build\n");
  return 1;
#else
  if (omp_get_max_threads() < 2) {
    scs_printf("skipping: OpenMP is limited to one thread\n");
    return 1;
  }
  return 0;
#endif
}

#endif /* PSD_PAR_BUILD */

#ifdef _OPENMP
/* largest difference between two projections that is rounding error */
#ifdef SFLOAT
#define PSD_PAR_TOL (1e-3)
#else
#define PSD_PAR_TOL (1e-9)
#endif

/* whether a and b agree: bit for bit if exact, else to rounding error */
static scs_int psd_par_same(const scs_float *a, const scs_float *b, scs_int m,
                            scs_int exact) {
  if (exact)
    return memcmp(a, b, m * sizeof(scs_float)) == 0;
  return SCS(norm_inf_diff)(a, b, m) <= PSD_PAR_TOL;
}

/* out = projection of v onto the dual cone in the metric r_y, with workspace c
 * and nt OpenMP threads. Returns the status of the projection. */
static scs_int psd_par_proj(ScsConeWork *c, const scs_float *v, scs_float *out,
                            scs_float *r_y, scs_int m, int nt) {
  memcpy(out, v, m * sizeof(scs_float));
  omp_set_num_threads(nt);
  return SCS(proj_dual_cone)(out, c, SCS_NULL, r_y);
}

/* 40 real blocks (four above LAPACK's blocking crossover, ~32, of different
 * sizes) and 4 complex ones: the workspaces are sized per thread, and if the
 * LAPACK calls were passed the lengths of their thread's work arrays instead
 * of the sizes of their own block, the large blocks would round differently,
 * and this mix shows it. The largest block holds less than half of the
 * work, so the parallel path is used. Other cones before and after the PSD
 * blocks check the offsets. */
static scs_int psd_par_s[] = {1,  2,  3, 5,  8, 13, 21, 1, 2,  3,  48, 5, 8, 13,
                              21, 1,  2, 3,  5, 8,  64, 21, 1, 2,  3,  5, 8, 13,
                              21, 1, 40, 3,  5, 8,  13, 36, 1, 2,  3,  5};
static scs_int psd_par_cs[] = {2, 3, 6, 10};
static scs_int psd_par_q[] = {4, 6};

/* number of rows of the cone k */
static scs_int psd_par_rows(const ScsCone *k) {
  scs_int i, m = k->z + k->l + 3 * k->ep;
  for (i = 0; i < k->qsize; ++i)
    m += k->q[i];
  for (i = 0; i < k->ssize; ++i)
    m += (k->s[i] * (k->s[i] + 1)) / 2;
  for (i = 0; i < k->cssize; ++i)
    m += k->cs[i] * k->cs[i];
  return m;
}

/* the cone of test_psd_parallel: its rows are returned */
static scs_int psd_par_cone(ScsCone *k) {
  memset(k, 0, sizeof(ScsCone));
  k->z = 5;
  k->l = 10;
  k->qsize = sizeof(psd_par_q) / sizeof(scs_int);
  k->q = psd_par_q;
  k->ssize = sizeof(psd_par_s) / sizeof(scs_int);
  k->s = psd_par_s;
  k->cssize = sizeof(psd_par_cs) / sizeof(scs_int);
  k->cs = psd_par_cs;
  k->ep = 2;
  return psd_par_rows(k);
}
#endif /* _OPENMP */

#ifdef PSD_PAR_BUILD
/* the cone of test_psd_parallel_solve: enough work in the PSD blocks for the
 * parallel path, also with an OpenBLAS whose lock makes blocks of n = 12 not
 * worth the threads, and no block with more than half of it */
static scs_int psd_par_solve_s[] = {24, 24, 22, 24};
static scs_int psd_par_solve_cs[] = {12};
static scs_int psd_par_solve_q[] = {5};

/* How close each solve's objectives, c'x and -b'y, must be to those of the
 * known solution: SOL_TOL, 1e-4 in double precision and 1e-2 in single.
 *
 * Single precision had the double-precision 1e-4, which a correct solve does
 * not meet. The harness raises the solves' eps 1e-5 to 1e-4
 * (problem_utils.h), and a solve stops once its gap is within eps_abs +
 * eps_rel |c'x| = 1.8e-3 here (the optimum is 17.0), so either objective can
 * be about that far off. Over 101 starting scales (moved by up to 50 ulps)
 * on the Mac, correct builds were off by up to 1.8e-3 with CG and 9e-4 with
 * QDLDL and the dense solver, and 23 to 67 of the starts failed 1e-4 in
 * every build (GCC 16 and clang, FMA contraction fast, on or off), among
 * them the test's own with CG in GCC 16's default build. The 4-thread and
 * 1-thread solves were the same bit for bit at every start, and so were
 * builds without the parallel projection or without OpenMP: the solve is
 * no less accurate in parallel. Double precision (eps 1e-5) was off by at
 * most 1.5e-5.
 *
 * What it catches: with the blocks of the last thread left unprojected, the
 * solve still reports solved, with an objective 4.2e-2 to 4.5e-2 off in
 * either precision (QDLDL, CG and dense), and verify_solution_correct passes,
 * since its cone distances go through the same projection (with the first
 * thread's blocks left unprojected it reports unbounded instead). 1e-2 is 4
 * times below that and 5.7 times above correct single-precision solves. */
#define PSD_PAR_OBJ_TOL SOL_TOL

/* Solve with nt OpenMP threads; psd_nt is set to the number of PSD threads of
 * the workspace. Returns NULL if the problem is solved and the solution is
 * correct. */
static const char *psd_par_solve(ScsData *d, ScsCone *k, ScsSettings *stgs,
                                 ScsSolution *opt_sol, ScsSolution *sol,
                                 ScsInfo *info, int nt, scs_int *psd_nt) {
  ScsWork *w;
  scs_int status;
  scs_float perr, derr;
  omp_set_num_threads(nt);
  w = scs_init(d, k, stgs);
  mu_assert("psd_par_solve: scs_init failed", w);
  *psd_nt = w->cone_work->psd_nthreads;
  status = scs_solve(w, sol, info, 0);
  scs_finish(w);
  perr = SCS(dot)(d->c, sol->x, d->n) - SCS(dot)(d->c, opt_sol->x, d->n);
  derr = -SCS(dot)(d->b, sol->y, d->m) + SCS(dot)(d->b, opt_sol->y, d->m);
  scs_printf("psd threads %li, primal error %4e, dual error %4e\n",
             (long)*psd_nt, perr, derr);
  mu_assert("psd_par_solve: SCS failed to produce SCS_SOLVED",
            status == SCS_SOLVED);
  mu_assert("psd_par_solve: objective off the known optimum",
            ABS(perr) < PSD_PAR_OBJ_TOL && ABS(derr) < PSD_PAR_OBJ_TOL);
  return verify_solution_correct(d, k, stgs, info, sol, status);
}
#endif /* PSD_PAR_BUILD */

static const char *test_psd_parallel(void) {
#ifndef _OPENMP
  return mu_skipped; /* the projection is serial without OpenMP */
#elif !defined(PSD_PAR_BUILD)
  /* OpenMP without USE_PARALLEL_PSD: serial, whatever the OpenMP threads */
  ScsCone k;
  scs_int m = psd_par_cone(&k), nt = 0;
  int omp_nt = omp_get_max_threads();
  ScsConeWork *c;
  omp_set_num_threads(4);
  c = SCS(init_cone)(&k, m);
  if (c)
    nt = c->psd_nthreads;
  SCS(finish_cone)(c);
  omp_set_num_threads(omp_nt);
  mu_assert("test_psd_parallel: cone setup failed", c);
  mu_assert("test_psd_parallel: PSD projection in parallel without "
            "USE_PARALLEL_PSD",
            nt == 1);
  scs_printf("psd blocks projected serially: no USE_PARALLEL_PSD\n");
  return 0;
#else
  ScsCone k;
  ScsConeWork *cs = SCS_NULL, *cp = SCS_NULL, *cc;
  scs_float *v, *vs, *vp, *r_y;
  scs_int i, m;
  const char *fail = SCS_NULL;
  PsdParState st;

  psd_par_enter(&st);
  if (psd_par_unavailable()) {
    psd_par_leave(&st);
    return mu_skipped;
  }

  m = psd_par_cone(&k);

  v = (scs_float *)scs_calloc(m, sizeof(scs_float));
  vs = (scs_float *)scs_calloc(m, sizeof(scs_float));
  vp = (scs_float *)scs_calloc(m, sizeof(scs_float));
  r_y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  if (v && vs && vp && r_y) {
    ran_start(99);
    for (i = 0; i < m; ++i) {
      v[i] = rand_scs_float();
      r_y[i] = 1.5 + rand_scs_float(); /* non-uniform metric */
    }
    omp_set_num_threads(1); /* an OpenMP BLAS then runs serially too */
    cs = SCS(init_cone)(&k, m);
    omp_set_num_threads(4);
    cp = SCS(init_cone)(&k, m);
  }

  if (!cs || !cp) {
    fail = "test_psd_parallel: cone setup failed";
  } else if (cs->psd_nthreads != 1 || cp->psd_nthreads < 2) {
    fail = "test_psd_parallel: parallel PSD projection not set up";
  } else if (psd_par_proj(cs, v, vs, r_y, m, 1) ||
             psd_par_proj(cp, v, vp, r_y, m, 4)) {
    fail = "test_psd_parallel: projection failed";
  } else if (!psd_par_same(vs, vp, m, st.exact)) {
    fail = "test_psd_parallel: parallel projection differs from serial";
  } else if (psd_par_proj(cp, v, vp, r_y, m, 2) ||
             !psd_par_same(vs, vp, m, st.exact)) {
    /* a team smaller than the split: each thread projects several lists */
    fail = "test_psd_parallel: 2 threads on a split for more differ";
  } else if (psd_par_proj(cp, v, vp, r_y, m, 1) ||
             !psd_par_same(vs, vp, m, st.exact)) {
    fail = "test_psd_parallel: 1 thread on a split for more differs";
  } else {
    scs_printf("psd blocks projected by %li threads, compared %s\n",
               (long)cp->psd_nthreads,
               st.exact ? "bit for bit" : "to rounding error");
  }

  if (!fail) {
    /* the OpenMP thread count at setup is the upper bound, also later */
    omp_set_num_threads(2);
    cc = SCS(init_cone)(&k, m);
    if (!cc || cc->psd_nthreads != 2)
      fail = "test_psd_parallel: 2 OpenMP threads at setup are not 2 threads";
    else if (psd_par_proj(cc, v, vp, r_y, m, 4) ||
             !psd_par_same(vs, vp, m, st.exact))
      fail = "test_psd_parallel: a split for 2 threads differs from serial";
    SCS(finish_cone)(cc);
  }

  SCS(finish_cone)(cs);
  SCS(finish_cone)(cp);
  scs_free(v);
  scs_free(vs);
  scs_free(vp);
  scs_free(r_y);
  psd_par_leave(&st);
  return fail;
#endif
}

/* No other solve in the tests has enough PSD blocks for the parallel path.
 * A random feasible problem over several cones, solved with 4 threads and
 * with 1: both must solve it, to the same iterates where the BLAS allows the
 * comparison. It has enough variables for every PSD block to matter: with
 * fewer (n = 50) a solver that left one list of blocks unprojected still
 * reached the optimum, and here it reports an unbounded problem. */
static const char *test_psd_parallel_solve(void) {
#ifndef PSD_PAR_BUILD
  return mu_skipped; /* the projection is serial in this build */
#else
  ScsCone k = {0};
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol1, *sol4, *opt_sol;
  ScsInfo info1 = {0}, info4 = {0};
  scs_int n = 200, col_nnz = 8, nt1 = 0, nt4 = 0;
  const char *fail;
  PsdParState st;

  psd_par_enter(&st);
  if (psd_par_unavailable()) {
    psd_par_leave(&st);
    return mu_skipped;
  }

  k.z = 4;
  k.l = 6;
  k.qsize = sizeof(psd_par_solve_q) / sizeof(scs_int);
  k.q = psd_par_solve_q;
  k.ssize = sizeof(psd_par_solve_s) / sizeof(scs_int);
  k.s = psd_par_solve_s;
  k.cssize = sizeof(psd_par_solve_cs) / sizeof(scs_int);
  k.cs = psd_par_solve_cs;
  k.ep = 2;

  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol1 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  sol4 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  d->m = psd_par_rows(&k);
  d->n = n;
  omp_set_num_threads(1);
  gen_random_prob_data(n * col_nnz, col_nnz, d, &k, opt_sol, 5678);
  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-5;
  stgs->eps_rel = 1e-5;

  fail = psd_par_solve(d, &k, stgs, opt_sol, sol4, &info4, 4, &nt4);
  if (!fail)
    fail = psd_par_solve(d, &k, stgs, opt_sol, sol1, &info1, 1, &nt1);
  if (!fail) {
    if (nt1 != 1 || nt4 < 2)
      fail = "test_psd_parallel_solve: parallel PSD projection not set up";
    else if (st.exact_solve &&
             (info1.iter != info4.iter ||
              !psd_par_same(sol1->x, sol4->x, d->n, 1) ||
              !psd_par_same(sol1->y, sol4->y, d->m, 1) ||
              !psd_par_same(sol1->s, sol4->s, d->m, 1)))
      fail = "test_psd_parallel_solve: iterates depend on the thread count";
  }

  SCS(free_data)(d);
  SCS(free_sol)(sol1);
  SCS(free_sol)(sol4);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  psd_par_leave(&st);
  return fail;
#endif
}

/* A build with USE_PARALLEL_PSD that projects the blocks one after another
 * sets the OpenMP thread count, and MKL's thread-local count, to one for the
 * blocks of 2 <= n <= PSD_SERIAL_HOLD_N and back after them (cones.c), so
 * that a BLAS on the OpenMP runtime of SCS, and MKL, project them on one
 * thread. The tests cannot see the threads of the BLAS. They check what the
 * hold must do on every platform: which blocks it holds (psd_nhold), on both
 * sides of the size, that the projection matches the one on one thread, and
 * that both counts are what they were after a projection and after a solve.
 * The held blocks are compared bit for bit with the reference on one thread
 * where the tests can tell the BLAS's rounding from the thread count (an
 * OpenBLAS on ELF). With an OpenMP OpenBLAS, which follows the OpenMP thread
 * count, that comparison is the one check here that sees the hold reach the
 * BLAS: without the hold, the blocks of 384 and less run on 4 threads and
 * round differently (the CI step with an OpenMP OpenBLAS in build_openmp.yml
 * runs it, and it fails there with the hold taken out). For MKL no check
 * here sees the hold act; the CI step with test/mkl_psd_guard.c in
 * build_openmp_mkl.yml does, from MKL's own report of its threads. No test
 * here sees the OpenMP thread count go back before a larger block within a
 * projection, only at its end: that block shows it only in the threads of
 * its BLAS (for MKL, that CI step sees it). Every other build holds nothing,
 * and its projections must leave the thread count alone too. */
#ifdef _OPENMP
/* PSD_SERIAL_HOLD_N of cones.c, which the tests cannot see */
#ifndef PSD_SERIAL_HOLD_N
#define PSD_SERIAL_HOLD_N (384)
#endif

/* The cones of test_psd_serial_hold, all projected serially: a free block
 * (above the size) first, then held blocks and a 1 x 1 block last, so that
 * the thread count goes back in the loop; a held block at the size itself
 * last, so that it goes back at the end; a held complex block; a free block
 * alone, so that nothing is held and MKL's count must be left alone: a
 * restore would write its value from before a hold that never happened (0)
 * over the 3 that psd_hold_check gives it on 2 threads. */
static scs_int psd_hold_s1[] = {PSD_SERIAL_HOLD_N + 1, 10, 3, 2, 1};
static scs_int psd_hold_s2[] = {1, PSD_SERIAL_HOLD_N};
static scs_int psd_hold_s3[] = {3};
static scs_int psd_hold_cs3[] = {200};
static scs_int psd_hold_s4[] = {PSD_SERIAL_HOLD_N + 1};

/* Blocks a serial projection of k holds: those of 2 <= n <= PSD_SERIAL_HOLD_N
 * in a build with the parallel projection, none in any other */
static scs_int psd_hold_expected(const ScsCone *k) {
  scs_int nh = 0;
#ifdef PSD_PAR_BUILD
  scs_int i;
  for (i = 0; i < k->ssize; ++i)
    nh += k->s[i] >= 2 && k->s[i] <= PSD_SERIAL_HOLD_N;
  for (i = 0; i < k->cssize; ++i)
    nh += k->cs[i] >= 2 && k->cs[i] <= PSD_SERIAL_HOLD_N;
#else
  (void)k;
#endif
  return nh;
}

/* the PSD blocks s and cs, with the other cones of psd_par_cone around them:
 * its rows are returned */
static scs_int psd_hold_cone(ScsCone *k, scs_int *s, scs_int ns, scs_int *cs,
                             scs_int ncs) {
  memset(k, 0, sizeof(ScsCone));
  k->z = 5;
  k->l = 10;
  k->qsize = sizeof(psd_par_q) / sizeof(scs_int);
  k->q = psd_par_q;
  k->ssize = ns;
  k->s = s;
  k->cssize = ncs;
  k->cs = cs;
  k->ep = 2;
  return psd_par_rows(k);
}

/* MKL's thread-local thread count of the calling thread (0: MKL's global
 * setting, where the tests leave it). A serial projection that sets it to
 * one for the blocks it holds (PSD_SERIAL_HOLD_N) must leave it as it was:
 * psd_mkl_local_set gives it a value before a projection and
 * psd_mkl_local_reset reads it after. Without MKL both do nothing. */

#ifdef PSD_PAR_BUILD
/* whether the tests can reach MKL's thread-local count */
static int psd_mkl_local_seen(void) {
#if PSD_PAR_WEAK_BLAS
  return MKL_Set_Num_Threads_Local ? 1 : 0;
#else
  return 0;
#endif
}
#endif

/* Set the count to nt, and return what psd_mkl_local_reset must then give:
 * nt with MKL, 0 without */
static int psd_mkl_local_set(int nt) {
#if defined(PSD_PAR_WEAK_BLAS) && PSD_PAR_WEAK_BLAS
  if (MKL_Set_Num_Threads_Local) {
    MKL_Set_Num_Threads_Local(nt);
    return nt;
  }
#endif
  (void)nt;
  return 0;
}

/* The count, which is then set back to 0; 0 without MKL */
static int psd_mkl_local_reset(void) {
#if defined(PSD_PAR_WEAK_BLAS) && PSD_PAR_WEAK_BLAS
  if (MKL_Set_Num_Threads_Local)
    return MKL_Set_Num_Threads_Local(0);
#endif
  return 0;
}

/* Project a random vector onto the dual cone of k on one OpenMP thread (the
 * reference) and with a cone set up on 4, on 4, 2 and 1 threads: the results
 * must agree (bit for bit if exact), and the OpenMP thread count and MKL's
 * thread-local count must be the same after each projection. MKL's starts
 * at 0 (its global setting) on 4 and 1 threads and at 3 on 2, so that a
 * restore to 0 or to the OpenMP thread count shows. On one thread the
 * projection is serial even if the cone was set up for the parallel loop,
 * and holds MKL for the same blocks. With serial, the cone must be set up
 * for a serial projection; nhold is set to the blocks it holds. */
static const char *psd_hold_check(ScsCone *k, scs_int serial, scs_int exact,
                                  scs_int *nhold) {
  static const int nts[] = {4, 2, 1}, mkl_start[] = {0, 3, 0};
  scs_int m = psd_par_rows(k), i, st;
  int nt, mkl_want, mkl_got;
  ScsConeWork *c1 = SCS_NULL, *c4 = SCS_NULL;
  const char *fail = SCS_NULL;
  scs_float *v = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *v1 = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *vh = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *r_y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  *nhold = 0;
  if (v && v1 && vh && r_y) {
    ran_start(7);
    for (i = 0; i < m; ++i) {
      v[i] = rand_scs_float();
      r_y[i] = 1.5 + rand_scs_float(); /* non-uniform metric */
    }
    omp_set_num_threads(1);
    c1 = SCS(init_cone)(k, m);
    omp_set_num_threads(4);
    c4 = SCS(init_cone)(k, m);
  }
  if (!c1 || !c4) {
    fail = "test_psd_serial_hold: cone setup failed";
  } else if (serial && c4->psd_nthreads != 1) {
    fail = "test_psd_serial_hold: the blocks are not projected serially";
  } else if (c4->psd_nhold != psd_hold_expected(k)) {
    fail = "test_psd_serial_hold: the wrong blocks are held";
  } else if (psd_par_proj(c1, v, v1, r_y, m, 1)) {
    fail = "test_psd_serial_hold: projection on one thread failed";
  } else {
    *nhold = c4->psd_nthreads == 1 ? c4->psd_nhold : 0;
    for (i = 0; i < 3 && !fail; ++i) {
      nt = nts[i];
      mkl_want = psd_mkl_local_set(mkl_start[i]);
      st = psd_par_proj(c4, v, vh, r_y, m, nt);
      mkl_got = psd_mkl_local_reset();
      if (st)
        fail = "test_psd_serial_hold: projection failed";
      else if (omp_get_max_threads() != nt)
        fail = "test_psd_serial_hold: OpenMP thread count not restored";
      else if (mkl_got != mkl_want)
        fail = "test_psd_serial_hold: MKL thread-local count not restored";
      else if (!psd_par_same(v1, vh, m, exact))
        fail = "test_psd_serial_hold: projection differs from one thread";
    }
  }
  SCS(finish_cone)(c1);
  SCS(finish_cone)(c4);
  scs_free(v);
  scs_free(v1);
  scs_free(vh);
  scs_free(r_y);
  return fail;
}
#endif /* _OPENMP */

static const char *test_psd_serial_hold(void) {
#ifndef _OPENMP
  return mu_skipped; /* no OpenMP thread count to hold */
#else
  ScsCone k;
  scs_int nh, nhold = 0, exact = 0, exact_free = 0;
  const char *fail;
#ifdef PSD_PAR_BUILD
  PsdParState st;
  psd_par_enter(&st);
  if (psd_par_unavailable()) {
    psd_par_leave(&st);
    return mu_skipped;
  }
  /* the held blocks run on one thread, as the reference does; the free one
   * follows the OpenMP thread count in an OpenMP OpenBLAS */
  exact = st.exact;
  exact_free = st.exact_solve;
#else
  int omp_nt = omp_get_max_threads();
  omp_set_num_threads(4);
  if (omp_get_max_threads() < 2) {
    omp_set_num_threads(omp_nt);
    scs_printf("skipping: OpenMP is limited to one thread\n");
    return mu_skipped;
  }
#endif

  psd_hold_cone(&k, psd_hold_s1, 5, SCS_NULL, 0);
  fail = psd_hold_check(&k, 1, exact_free, &nh);
  nhold += nh;
  if (!fail) {
    psd_hold_cone(&k, psd_hold_s2, 2, SCS_NULL, 0);
    fail = psd_hold_check(&k, 1, exact, &nh);
    nhold += nh;
  }
  if (!fail) {
    psd_hold_cone(&k, psd_hold_s3, 1, psd_hold_cs3, 1);
    fail = psd_hold_check(&k, 1, exact, &nh);
    nhold += nh;
  }
  if (!fail) {
    psd_hold_cone(&k, psd_hold_s4, 1, SCS_NULL, 0);
    fail = psd_hold_check(&k, 1, exact_free, &nh);
    nhold += nh;
  }
  if (!fail) {
    /* the cone of test_psd_parallel: set up for the parallel loop where the
     * build has it, and then its blocks are held only on one thread */
    psd_par_cone(&k);
    fail = psd_hold_check(&k, 0, exact, &nh);
  }

#ifdef PSD_PAR_BUILD
  psd_par_leave(&st);
#else
  omp_set_num_threads(omp_nt);
#endif
  if (fail)
    return fail;
#ifdef PSD_PAR_BUILD
  mu_assert("test_psd_serial_hold: no block held", nhold > 0);
  scs_printf("psd serial hold: %li blocks of 2 <= n <= %li held, OpenMP "
             "thread count%s restored, compared %s\n",
             (long)nhold, (long)PSD_SERIAL_HOLD_N,
             psd_mkl_local_seen() ? " and MKL thread-local count" : "",
             exact ? "bit for bit" : "to rounding error");
#else
  mu_assert("test_psd_serial_hold: blocks held without USE_PARALLEL_PSD",
            nhold == 0);
  scs_printf("psd serial hold: none, no USE_PARALLEL_PSD\n");
#endif
  return 0;
#endif
}

/* A solve whose blocks are projected serially, some of them held, on 4 and
 * on 1 OpenMP threads: both must solve it, to the same iterates where the
 * BLAS allows the comparison, and leave the OpenMP thread count as it was. */
static const char *test_psd_serial_hold_solve(void) {
#ifndef PSD_PAR_BUILD
  return mu_skipped; /* nothing is held in this build */
#else
  static scs_int s[] = {30, 3, 2};
  ScsCone k = {0};
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol1, *sol4, *opt_sol;
  ScsInfo info1 = {0}, info4 = {0};
  scs_int n = 100, col_nnz = 8, nt1 = 0, nt4 = 0;
  int mkl_want, mkl_got;
  const char *fail;
  PsdParState st;

  psd_par_enter(&st);
  if (psd_par_unavailable()) {
    psd_par_leave(&st);
    return mu_skipped;
  }

  k.z = 4;
  k.l = 6;
  k.qsize = sizeof(psd_par_solve_q) / sizeof(scs_int);
  k.q = psd_par_solve_q;
  k.ssize = sizeof(s) / sizeof(scs_int);
  k.s = s;
  k.ep = 2;

  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol1 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  sol4 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  d->m = psd_par_rows(&k);
  d->n = n;
  omp_set_num_threads(1);
  gen_random_prob_data(n * col_nnz, col_nnz, d, &k, opt_sol, 4321);
  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-5;
  stgs->eps_rel = 1e-5;

  /* MKL's thread-local count starts at 3 for the solve on 4 threads and at 0
   * (MKL's global setting) for the one on 1 thread */
  mkl_want = psd_mkl_local_set(3);
  fail = psd_par_solve(d, &k, stgs, opt_sol, sol4, &info4, 4, &nt4);
  mkl_got = psd_mkl_local_reset();
  if (!fail && omp_get_max_threads() != 4)
    fail = "test_psd_serial_hold_solve: OpenMP thread count not restored";
  if (!fail && mkl_got != mkl_want)
    fail = "test_psd_serial_hold_solve: MKL thread-local count not restored";
  if (!fail) {
    fail = psd_par_solve(d, &k, stgs, opt_sol, sol1, &info1, 1, &nt1);
    mkl_got = psd_mkl_local_reset();
  }
  if (!fail && omp_get_max_threads() != 1)
    fail = "test_psd_serial_hold_solve: OpenMP thread count not restored";
  if (!fail && mkl_got != 0)
    fail = "test_psd_serial_hold_solve: MKL thread-local count not restored";
  if (!fail) {
    if (nt1 != 1 || nt4 != 1)
      fail = "test_psd_serial_hold_solve: the blocks are not projected "
             "serially";
    else if (st.exact_solve &&
             (info1.iter != info4.iter ||
              !psd_par_same(sol1->x, sol4->x, d->n, 1) ||
              !psd_par_same(sol1->y, sol4->y, d->m, 1) ||
              !psd_par_same(sol1->s, sol4->s, d->m, 1)))
      fail = "test_psd_serial_hold_solve: iterates depend on the thread "
             "count";
  }

  SCS(free_data)(d);
  SCS(free_sol)(sol1);
  SCS(free_sol)(sol4);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);
  psd_par_leave(&st);
  return fail;
#endif
}
