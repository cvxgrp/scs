/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */
#include <stdio.h>

#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"

/* Include Tests */
#include "problems/rob_gauss_cov_est.h"
#include "problems/small_lp.h"

int tests_run = 0;

static const char *all_tests(void) {
  mu_run_test(small_lp);
  mu_run_test(rob_gauss_cov_est);
  return 0;
}

int main(void) {
  const char *result = all_tests();
  if (result != 0) {
    printf("%s\n", result);
  } else {
    printf("ALL TESTS PASSED\n");
  }
  printf("Tests run: %d\n", tests_run);

  return result != 0;
}
