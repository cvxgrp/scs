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

/* With OpenMP the PSD blocks are projected in parallel (cones.c). Each block
 * goes through the same LAPACK calls with the same workspace sizes as in the
 * serial projection, so with single-threaded BLAS the result is the same bit
 * for bit, whatever the number of threads. test_psd_parallel projects one
 * vector with 1 and with 4 threads, and test_psd_parallel_solve solves one
 * problem both ways.
 *
 * The comparisons are bit for bit when the BLAS is an OpenBLAS that the tests
 * can set to one thread. They are to rounding error for the other libraries,
 * which may round one block differently from call to call (MKL: the result
 * can depend on the alignment of the workspace) or cannot be set to one
 * thread from here (Accelerate). */

#ifdef _OPENMP
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
#ifdef __cplusplus
}
#endif
#else
#define PSD_PAR_WEAK_BLAS (0) /* e.g. Accelerate: VECLIB_MAXIMUM_THREADS=1 */
#endif

/* largest difference between two projections that is rounding error */
#ifdef SFLOAT
#define PSD_PAR_TOL (1e-3)
#else
#define PSD_PAR_TOL (1e-9)
#endif

#define PSD_PAR_ENV "SCS_PSD_NUM_THREADS"

/* the thread settings before a test, and what the test can rely on */
typedef struct {
  int omp_nt;          /* OpenMP threads */
  int ob_nt;           /* OpenBLAS threads, 0: no OpenBLAS found */
  int mkl_nt;          /* MKL threads, 0: no MKL found */
  char *env;           /* SCS_PSD_NUM_THREADS, NULL: not set */
  scs_int exact;       /* projections are the same bit for bit */
  scs_int exact_solve; /* and so is everything else in a solve */
} PsdParState;

/* set SCS_PSD_NUM_THREADS to val, or unset it if val is NULL */
static void psd_par_setenv(const char *val) {
#ifdef _WIN32
  _putenv_s(PSD_PAR_ENV, val ? val : "");
#else
  if (val)
    setenv(PSD_PAR_ENV, val, 1);
  else
    unsetenv(PSD_PAR_ENV);
#endif
}

/* Save the thread settings in st. Then no SCS_PSD_NUM_THREADS, and BLAS on one
 * thread where the library can be told: a threaded BLAS rounds differently in
 * the serial projection, and a pthreads OpenBLAS with more than one thread
 * keeps the projection serial. */
static void psd_par_enter(PsdParState *st) {
  const char *env = getenv(PSD_PAR_ENV);
  memset(st, 0, sizeof(PsdParState));
  st->omp_nt = omp_get_max_threads();
  if (env) {
    st->env = (char *)scs_malloc(strlen(env) + 1);
    if (st->env)
      memcpy(st->env, env, strlen(env) + 1);
  }
  psd_par_setenv(SCS_NULL);
#if PSD_PAR_WEAK_BLAS
  if (openblas_set_num_threads && openblas_get_num_threads) {
    st->ob_nt = openblas_get_num_threads();
    openblas_set_num_threads(1);
    st->exact = !MKL_Set_Num_Threads;
    /* an OpenMP OpenBLAS (2) follows the OpenMP thread count in the BLAS calls
     * outside of the PSD threads */
    st->exact_solve =
        st->exact && openblas_get_parallel && openblas_get_parallel() != 2;
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
  psd_par_setenv(st->env);
  scs_free(st->env);
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
 * sizes) and 4 complex ones: with workspaces sized per thread instead of the
 * serial sizes, the large blocks would round differently, and this mix shows
 * it. The largest block holds less than half of the work, so the parallel
 * path is used. Other cones before and after the PSD blocks check the
 * offsets. */
static scs_int psd_par_s[] = {1,  2,  3, 5,  8, 13, 21, 1, 2,  3,  48, 5, 8, 13,
                              21, 1,  2, 3,  5, 8,  64, 21, 1, 2,  3,  5, 8, 13,
                              21, 1, 40, 3,  5, 8,  13, 36, 1, 2,  3,  5};
static scs_int psd_par_cs[] = {2, 3, 6, 10};
static scs_int psd_par_q[] = {4, 6};

/* the cone of test_psd_parallel_solve: enough work in the PSD blocks for the
 * parallel path, and no block with more than half of it */
static scs_int psd_par_solve_s[] = {12, 12, 11, 12};
static scs_int psd_par_solve_cs[] = {6};
static scs_int psd_par_solve_q[] = {5};

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
  mu_assert("test_psd_parallel_solve: scs_init failed", w);
  *psd_nt = w->cone_work->psd_nthreads;
  status = scs_solve(w, sol, info, 0);
  scs_finish(w);
  perr = SCS(dot)(d->c, sol->x, d->n) - SCS(dot)(d->c, opt_sol->x, d->n);
  derr = -SCS(dot)(d->b, sol->y, d->m) + SCS(dot)(d->b, opt_sol->y, d->m);
  scs_printf("psd threads %li, primal error %4e, dual error %4e\n",
             (long)*psd_nt, perr, derr);
  mu_assert("test_psd_parallel_solve: SCS failed to produce SCS_SOLVED",
            status == SCS_SOLVED && ABS(perr) < 1e-4 && ABS(derr) < 1e-4);
  return verify_solution_correct(d, k, stgs, info, sol, status);
}
#endif /* _OPENMP */

static const char *test_psd_parallel(void) {
#ifndef _OPENMP
  return mu_skipped; /* the projection is serial without OpenMP */
#else
  ScsCone k = {0};
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

  k.z = 5;
  k.l = 10;
  k.qsize = sizeof(psd_par_q) / sizeof(scs_int);
  k.q = psd_par_q;
  k.ssize = sizeof(psd_par_s) / sizeof(scs_int);
  k.s = psd_par_s;
  k.cssize = sizeof(psd_par_cs) / sizeof(scs_int);
  k.cs = psd_par_cs;
  k.ep = 2;
  m = psd_par_rows(&k);

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
    /* SCS_PSD_NUM_THREADS, read at setup, is an upper bound */
    omp_set_num_threads(4);
    psd_par_setenv("1");
    cc = SCS(init_cone)(&k, m);
    if (!cc || cc->psd_nthreads != 1)
      fail = "test_psd_parallel: SCS_PSD_NUM_THREADS=1 is not serial";
    SCS(finish_cone)(cc);
    psd_par_setenv("2");
    cc = SCS(init_cone)(&k, m);
    if (!fail && (!cc || cc->psd_nthreads != 2))
      fail = "test_psd_parallel: SCS_PSD_NUM_THREADS=2 is not 2 threads";
    else if (!fail && (psd_par_proj(cc, v, vp, r_y, m, 4) ||
                       !psd_par_same(vs, vp, m, st.exact)))
      fail = "test_psd_parallel: SCS_PSD_NUM_THREADS=2 differs from serial";
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
#ifndef _OPENMP
  return mu_skipped; /* the projection is serial without OpenMP */
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
