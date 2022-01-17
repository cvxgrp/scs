/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */
#include <stdio.h>

#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"

/* Include Tests */
#include "problems/degenerate.h"
#include "problems/hs21_tiny_qp.h"
#include "problems/hs21_tiny_qp_rw.h"
#include "problems/infeasible_tiny_qp.h"
#include "problems/qafiro_tiny_qp.h"
#include "problems/small_lp.h"
#include "problems/small_qp.h"
#include "problems/test_validation.h"
#include "problems/unbounded_tiny_qp.h"

#define _SKIP(problem)                                                         \
  char *problem(void) {                                                        \
    scs_printf("skipped\n");                                                   \
    return 0;                                                                  \
  }

#ifdef USE_LAPACK /* solve SDPs, requires blas / lapack */
#include "problems/rob_gauss_cov_est.h"
#else
_SKIP(rob_gauss_cov_est)
#endif

/* TODO: this reads a file written with 32bit ints */
#if defined(USE_LAPACK) && !defined(DLONG)
#include "problems/random_prob.h"
#else
_SKIP(random_prob)
#endif

int tests_run = 0;

static const char *all_tests(void) {
  mu_run_test(test_validation);
  mu_run_test(degenerate);
  mu_run_test(small_lp);
  mu_run_test(small_qp);
  mu_run_test(rob_gauss_cov_est);
  mu_run_test(hs21_tiny_qp);
  mu_run_test(hs21_tiny_qp_rw);
  mu_run_test(qafiro_tiny_qp);
  mu_run_test(infeasible_tiny_qp);
  mu_run_test(unbounded_tiny_qp);
  mu_run_test(random_prob);
  return 0;
}

int main(void) {
  const char *result = all_tests();
  if (result != 0) {
    scs_printf("%s\n", result);
  } else {
    scs_printf("ALL TESTS PASSED\n");
  }
  scs_printf("Tests run: %d\n", tests_run);

  return result != 0;
}
