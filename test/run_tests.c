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
#include "problems/rob_gauss_cov_est.h"
#include "problems/small_lp.h"
#include "problems/unbounded_tiny_qp.h"

int tests_run = 0;

static const char *all_tests(void) {
  scs_printf("degenerate\n");
  mu_run_test(degenerate);
  scs_printf("small_lp\n");
  mu_run_test(small_lp);
  scs_printf("rob_gauss_cov_est\n");
  mu_run_test(rob_gauss_cov_est);
  scs_printf("hs21_tiny_qp\n");
  mu_run_test(hs21_tiny_qp);
  scs_printf("hs21_tiny_qp_rw\n");
  mu_run_test(hs21_tiny_qp_rw);
  scs_printf("qafiro_tiny_qp\n");
  mu_run_test(qafiro_tiny_qp);
  scs_printf("infeasible_tiny_qp\n");
  mu_run_test(infeasible_tiny_qp);
  scs_printf("unbounded_tiny_qp\n");
  mu_run_test(unbounded_tiny_qp);
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
