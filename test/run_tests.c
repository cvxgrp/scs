/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */
#include <stdio.h>

#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"

/* Include Tests */
#include "problems/degenerate.h"
#include "problems/hs21_tiny_qp.h"
#include "problems/infeasible_tiny_qp.h"
#include "problems/qafiro_tiny_qp.h"
#include "problems/small_lp.h"
#include "problems/small_qp.h"
#include "problems/unbounded_tiny_qp.h"

int tests_run = 0;

/* decrement tests_run since mu_unit will increment it, so this cancels */
#define _SKIP(problem)                                                         \
  char *problem(void) {                                                        \
    scs_printf("skipped\n");                                                   \
    tests_run--;                                                               \
    return 0;                                                                  \
  }

#if NO_VALIDATE == 0
#include "problems/test_validation.h"
#else
_SKIP(test_validation)
#endif

/* solve SDPs, requires blas / lapack */
#if defined(USE_LAPACK) && NO_READ_WRITE == 0
#include "problems/rob_gauss_cov_est.h"
#else
_SKIP(rob_gauss_cov_est)
#endif

#if NO_READ_WRITE == 0 /* reads / writes */
#include "problems/hs21_tiny_qp_rw.h"
#else
_SKIP(hs21_tiny_qp_rw)
#endif

/* TODO: this reads a file written with 32bit ints */
#if defined(USE_LAPACK) && DLONG == 0 && NO_READ_WRITE == 0
#include "problems/random_prob.h"
#else
_SKIP(random_prob)
#endif

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
