#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

#include <math.h>
#include <string.h>

/*
 * Problems from public benchmark collections, solved as the solver_benchmarks
 * harness poses them for SCS. Each covers a case the synthetic tests do not,
 * on real data (see each problem below). The data are in problems/benchmark/,
 * one generated header per problem; README.md there gives each problem's
 * source, how it was posed and its licence, which SCS does not control.
 *
 * Each test solves its problem once from SCS's default settings, with
 * eps_abs = eps_rel set per problem, the default eps_infeas, and a max_iters
 * far above what it takes (at least 13 times the most iterations seen in the
 * screen below), and checks:
 *   - the exit flag, exactly: no inaccurate status passes;
 *   - verify_solution_correct (problem_utils.h): SCS's residuals,
 *     objectives and certificate conditions, recomputed from the data;
 *   - a solved problem: its primal and dual objectives against the
 *     reference optimum, to 1e-3 (1 + |opt|) in double precision and 1e-2
 *     in single;
 *   - the QPs: the original bounds, checked independently of SCS's
 *     residuals (bench_bounds_violation).
 * The tolerances are set by what a bug would give, not by the accuracy seen
 * locally. Iteration counts are printed, never checked: they change with the
 * backend, the BLAS, the compiler and eigCG.
 *
 * In single precision the data round to float, and the suite raises eps and
 * eps_infeas to 1e-4 (problem_utils.h). A problem is skipped only where the
 * screen found a failure or a cost that rules it out; the reason is printed.
 * The SDPs need USE_LAPACK for the PSD cone (run_tests.c skips them
 * without it).
 *
 * Screen (Mac, on master of 2026-11-25): each test's own settings and
 * checks from the default scale and 40 one-ulp moves of it (200 for CG and
 * single precision), on every backend in 24 builds (clang and gcc, -O0 to
 * -O3, Accelerate and OpenBLAS, 4- and 8-byte integers, with and without
 * LAPACK, eigCG on and off, ASan and UBSan): 59,009 solves, none failed.
 * Each problem's comment gives the iterations seen.
 */

/* skip flags: the builds a problem does not run in */
#define BENCH_SKIP_SFLOAT (1)          /* single precision */
#define BENCH_SKIP_INDIRECT (2)        /* the CG solvers, CPU and GPU */
#define BENCH_SKIP_SFLOAT_INDIRECT (4) /* single precision with CG */

typedef struct {
  const char *name;                    /* as its collection names it */
  void (*set)(ScsData *d, ScsCone *k); /* from its generated header */
  scs_int exitflag;                    /* expected, exactly */
  scs_float opt;                       /* optimal value, if SCS_SOLVED */
  scs_float eps;                       /* eps_abs = eps_rel */
  scs_int max_iters;
  /* > 0: check the original bounds (bench_bounds_violation) to this in
   * double precision, and to 1e-2 in single */
  scs_float bounds_tol;
  int skip; /* BENCH_SKIP_* */
  const char *skip_why;
} BenchProblem;

/* What one solve gave */
typedef struct {
  scs_int exitflag;
  scs_int iter;
  scs_float pobj_err, dobj_err; /* |objective - opt| / (1 + |opt|) */
  scs_float obj_tol;            /* 0 if not checked */
  scs_float bounds;             /* bench_bounds_violation */
  scs_float bounds_tol;         /* 0 if not checked */
} BenchResult;

/*
 * How far x is from meeting the original rows, with the box cone's t fixed
 * to 1: the largest violation over the rows, divided by 1 + the largest of
 * |a'x| and the finite |b|, |bl|, |bu|. This is the harness's relative
 * primal residual of l <= A x <= u (zero-cone rows are its equalities,
 * box rows its inequalities); the rows it leaves out, which have no finite
 * bound, would only add to the denominator.
 *
 * Why: the box cone poses l <= a'x <= u as (t, s) with t l <= s <= t u,
 * s = a'x, and a first row that says t = 1. SCS's own residuals, and
 * verify_solution_correct, measure that form, where t is a variable, so they
 * pass a point with t off 1 whose rows meet their bounds only up to |t - 1|
 * times the bound. SCS 3.3.2, which weighted the t row like an inequality,
 * returned such points as solved on GOULDQP2. Since #548 the t row is
 * weighted like an equality (test_box_cone_t_weight.h). Computed in double
 * whatever the precision of the build.
 */
static double bench_bounds_violation(const ScsData *d, const ScsCone *k,
                                     const scs_float *x) {
  scs_int i, j, p, box0 = k->z + k->l + 1; /* the first box row after t */
  double *ax = (double *)scs_calloc(d->m, sizeof(double));
  double viol = 0., size = 0., v, b, bl, bu;
  for (j = 0; j < d->n; ++j) {
    for (p = d->A->p[j]; p < d->A->p[j + 1]; ++p) {
      ax[d->A->i[p]] += (double)d->A->x[p] * (double)x[j];
    }
  }
  /* zero cone: A x = b; nonnegative cone: A x <= b */
  for (i = 0; i < k->z + k->l; ++i) {
    b = (double)d->b[i];
    viol = MAX(viol, i < k->z ? fabs(ax[i] - b) : MAX(ax[i] - b, 0.));
    size = MAX(size, MAX(fabs(ax[i]), fabs(b)));
  }
  /* box cone after t: bl <= b - A x <= bu, where b - A x = a'x */
  for (j = 0; j < k->bsize - 1; ++j) {
    i = box0 + j;
    v = (double)d->b[i] - ax[i];
    bl = (double)k->bl[j];
    bu = (double)k->bu[j];
    size = MAX(size, fabs(v));
    if (isfinite(bl)) {
      size = MAX(size, fabs(bl));
    }
    if (isfinite(bu)) {
      size = MAX(size, fabs(bu));
    }
    viol = MAX(viol, (isfinite(bl) ? MAX(bl - v, 0.) : 0.) +
                         (isfinite(bu) ? MAX(v - bu, 0.) : 0.));
  }
  scs_free(ax);
  return viol / (1. + size);
}

/* Whether this build skips b */
static int bench_skipped(const BenchProblem *b) {
  int sfloat = 0;
  int indirect = strstr(scs_get_lin_sys_method(), "indirect") != SCS_NULL;
#ifdef SFLOAT
  sfloat = 1;
#endif
  return ((b->skip & BENCH_SKIP_SFLOAT) && sfloat) ||
         ((b->skip & BENCH_SKIP_INDIRECT) && indirect) ||
         ((b->skip & BENCH_SKIP_SFLOAT_INDIRECT) && sfloat && indirect);
}

/* The settings each test solves with */
static void bench_settings(const BenchProblem *b, ScsSettings *stgs) {
  scs_set_default_settings(stgs);
  stgs->eps_abs = b->eps;
  stgs->eps_rel = b->eps;
  stgs->max_iters = b->max_iters;
  stgs->verbose = 1;
}

/* The checks of a solved problem, after verify_solution_correct */
static const char *bench_check_solved(const BenchProblem *b, ScsData *d,
                                      ScsCone *k, ScsInfo *info,
                                      ScsSolution *sol, BenchResult *r) {
  r->obj_tol = VERIFY_TOL(1e-3, 1e-2);
  r->pobj_err = ABS(info->pobj - b->opt) / (1 + ABS(b->opt));
  r->dobj_err = ABS(info->dobj - b->opt) / (1 + ABS(b->opt));
  mu_assert_less("primal objective off the optimum", r->pobj_err,
                 r->obj_tol);
  mu_assert_less("dual objective off the optimum", r->dobj_err, r->obj_tol);
  if (b->bounds_tol > 0) {
    r->bounds_tol = VERIFY_TOL(b->bounds_tol, 1e-2);
    r->bounds = (scs_float)bench_bounds_violation(d, k, sol->x);
    mu_assert_less("original bounds not met", r->bounds, r->bounds_tol);
  }
  return SCS_NULL;
}

/* Solve b with stgs and run the checks; r gets what the solve gave */
static const char *bench_solve(const BenchProblem *b, ScsSettings *stgs,
                               BenchResult *r) {
  ScsData d;
  ScsCone k;
  ScsSolution sol = {SCS_NULL, SCS_NULL, SCS_NULL};
  ScsInfo info = {0};
  const char *fail = SCS_NULL;

  memset(r, 0, sizeof(*r));
  b->set(&d, &k);
  r->exitflag = scs(&d, &k, stgs, &sol, &info);
  r->iter = info.iter;
  if (r->exitflag != b->exitflag) {
    scs_printf("%s: exit flag %li, expected %li\n", b->name,
               (long)r->exitflag, (long)b->exitflag);
    fail = "wrong exit flag";
  }
  if (!fail) {
    fail = verify_solution_correct(&d, &k, stgs, &info, &sol, r->exitflag);
  }
  if (!fail && b->exitflag == SCS_SOLVED) {
    fail = bench_check_solved(b, &d, &k, &info, &sol, r);
  }
  scs_printf("%s: %s, %li iterations", b->name, info.status, (long)r->iter);
  if (r->obj_tol > 0) {
    scs_printf(", objective off by %.1e (limit %.0e)",
               MAX(r->pobj_err, r->dobj_err), r->obj_tol);
  }
  if (r->bounds_tol > 0) {
    scs_printf(", bounds off by %.1e (limit %.0e)", r->bounds, r->bounds_tol);
  }
  scs_printf("\n");
  scs_free(sol.x);
  scs_free(sol.y);
  scs_free(sol.s);
  return fail;
}

static const char *bench_run(const BenchProblem *b) {
  ScsSettings stgs;
  BenchResult r;
  if (bench_skipped(b)) {
    scs_printf("%s: skipped: %s\n", b->name, b->skip_why);
    return mu_skipped;
  }
  bench_settings(b, &stgs);
  return bench_solve(b, &stgs, &r);
}

#define BENCH_TEST(id)                                                         \
  static const char *bench_##id(void) {                                        \
    return bench_run(&bench_##id##_problem);                                   \
  }

#if defined(USE_LAPACK)
#include "benchmark/infd2.h"
#include "benchmark/infp1.h"
#include "benchmark/qap5.h"
#include "benchmark/truss4.h"

/*
 * SDPLIB infp1: an SDP infeasibility certificate, one 30 x 30 PSD block
 * (m = 465, n = 10, A dense). No other test has an SDP certificate.
 * Screen: 50 iterations everywhere.
 */
static const BenchProblem bench_infp1_problem = {
    "infp1",
    bench_infp1_set,
    SCS_INFEASIBLE,
    0.,
    1e-4,
    1000,
    0.,
    0,
    SCS_NULL};

/*
 * SDPLIB infd2: an SDP unboundedness certificate, one 30 x 30 block. In
 * double precision it also turns on per-iteration refinement of the g solve.
 * Screen: 75 iterations in double precision, 75 to 300 in single.
 * Skipped in single precision with CG, where one of 201 starts did not
 * finish in 6,500 iterations.
 */
static const BenchProblem bench_infd2_problem = {
    "infd2",
    bench_infd2_set,
    SCS_UNBOUNDED,
    0.,
    1e-4,
    5000,
    0.,
    BENCH_SKIP_SFLOAT_INDIRECT,
    "single precision with CG can take thousands of iterations"};

/*
 * SDPLIB truss4: a solved SDP with many small blocks (six 3 x 3, one 1 x 1).
 * Optimal value -9.0099963 (SDPLIB). Screen: 125 iterations in double
 * precision, 125 to 150 in single.
 */
static const BenchProblem bench_truss4_problem = {
    "truss4",
    bench_truss4_set,
    SCS_SOLVED,
    -9.0099963,
    1e-4,
    2500,
    0.,
    0,
    SCS_NULL};

/*
 * SDPLIB qap5: a solved SDP with one 26 x 26 block, at eps 1e-6, which also
 * exercises CG and eigCG on a real SDP (eigCG changes its CG work).
 * Optimal value -436 (SDPLIB). Screen: in double precision, 275 iterations
 * on the direct solvers, and on CG 275 to 300 with eigCG and 325 to 550
 * without it (DEFLATE_VECTORS=0); 250 to 325 in single. Skipped in single
 * precision with CG, where about a third of 201 starts did not finish in
 * 6,750 iterations.
 */
static const BenchProblem bench_qap5_problem = {
    "qap5",
    bench_qap5_set,
    SCS_SOLVED,
    -436.,
    1e-6,
    20000,
    0.,
    BENCH_SKIP_SFLOAT_INDIRECT,
    "single precision with CG does not converge on every start"};

BENCH_TEST(infp1)
BENCH_TEST(infd2)
BENCH_TEST(truss4)
BENCH_TEST(qap5)
#endif

#include "benchmark/cvxqp2_s.h"
#include "benchmark/dualc1.h"
#include "benchmark/galenet.h"
#include "benchmark/gouldqp2.h"
#include "benchmark/itest2.h"

/*
 * Netlib lp/infeas galenet: an LP infeasibility certificate on real data,
 * through the box cone, with a zero cone (2 equality rows, 14 box rows). The
 * other LP certificate tests are toys.
 * Screen: 50 to 75 iterations.
 */
static const BenchProblem bench_galenet_problem = {
    "galenet",
    bench_galenet_set,
    SCS_INFEASIBLE,
    0.,
    1e-4,
    1000,
    0.,
    0,
    SCS_NULL};

/*
 * Netlib lp/infeas itest2: the same with no zero cone (all 13 rows in the box
 * cone), at eps 1e-6.
 * Screen: 50 to 75 iterations.
 */
static const BenchProblem bench_itest2_problem = {
    "itest2",
    bench_itest2_set,
    SCS_INFEASIBLE,
    0.,
    1e-6,
    1000,
    0.,
    0,
    SCS_NULL};

/*
 * Maros-Meszaros CVXQP2_S: a QP with P through the box cone at eps 1e-6,
 * where the t weight matters (3.3.2 took 225 iterations, 3.4.0 takes 75).
 * Optimal value 8120.9404773 (Clarabel at 1e-10; 00README.QP
 * 8.1209405e+03). Screen: 75 to 125 iterations.
 */
static const BenchProblem bench_cvxqp2_s_problem = {
    "CVXQP2_S",
    bench_cvxqp2_s_set,
    SCS_SOLVED,
    8120.9404773,
    1e-6,
    2500,
    1e-3,
    0,
    SCS_NULL};

/*
 * Maros-Meszaros DUALC1: CG and eigCG on a real QP (9 variables, 1 equality,
 * 223 box rows) whose CG work eigCG changes.
 * Optimal value 6155.2508295 (Clarabel at 1e-10; 00README.QP 6.1552508e+03).
 * Screen: 125 iterations. Double precision only: in single precision 12 of
 * 2,211 starts failed verify_solution_correct (the reported dual residual
 * against the recomputed one).
 */
static const BenchProblem bench_dualc1_problem = {
    "DUALC1",
    bench_dualc1_set,
    SCS_SOLVED,
    6155.2508295,
    1e-4,
    2500,
    1e-3,
    BENCH_SKIP_SFLOAT,
    "single precision fails verify_solution_correct on some starts"};

/*
 * Maros-Meszaros GOULDQP2: guards the box cone's accuracy, with the bounds
 * check at 5e-4. SCS 3.3.2 returns points as solved that pass
 * verify_solution_correct but miss the original bounds by up to 3.2e-3
 * (relative): this test fails 3.3.2 at 139 of 164 starts (QDLDL, dense and
 * Accelerate, and QDLDL without LAPACK), the default start among them on
 * QDLDL and Accelerate. The worst seen on 3.4.0 is 1.0e-4 (3,500 solves), 5
 * times inside the limit.
 * Optimal value 1.8427450e-4 (Clarabel at 1e-10; 00README.QP 1.8427534e-04).
 * Screen: 125 to 150 iterations. Double precision only, where the check can
 * tell 3.3.2's error from rounding, and not with CG, which takes 0.3 s here
 * (65,000 to 89,000 CG iterations).
 */
static const BenchProblem bench_gouldqp2_problem = {
    "GOULDQP2",
    bench_gouldqp2_set,
    SCS_SOLVED,
    1.8427450e-4,
    1e-4,
    2500,
    5e-4,
    BENCH_SKIP_SFLOAT | BENCH_SKIP_INDIRECT,
    "double precision only, and not with CG (0.3 s)"};

BENCH_TEST(galenet)
BENCH_TEST(itest2)
BENCH_TEST(cvxqp2_s)
BENCH_TEST(dualc1)
BENCH_TEST(gouldqp2)
