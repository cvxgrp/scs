/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */

/* Simple Macros for testing */

/* Failures are recorded rather than aborting the run, so that a single run
 * reports every failing test. The runner defines these. */
#define MU_MAX_FAILED_TESTS 256
extern int tests_run;
extern int tests_failed;
extern const char *failed_tests[MU_MAX_FAILED_TESTS];
#define mu_assert_less(message, a, b)                                          \
  do {                                                                         \
    if (a > b) {                                                               \
      scs_printf("%s: %1.3e > %1.3e\n", message, a, b);                        \
      return message;                                                          \
    }                                                                          \
  } while (0)

#define mu_assert(message, test)                                               \
  do {                                                                         \
    if (!(test))                                                               \
      return message;                                                          \
  } while (0)

#define mu_run_test(test) _mu_run_test(#test, test)

#define _mu_run_test(name, test)                                               \
  do {                                                                         \
    scs_printf("*********************************************************\n"); \
    scs_printf("Running test: %s\n", name);                                    \
    const char *message = test();                                              \
    tests_run++;                                                               \
    if (message) {                                                             \
      scs_printf("FAILED: %s: %s\n", name, message);                           \
      if (tests_failed < MU_MAX_FAILED_TESTS)                                  \
        failed_tests[tests_failed] = name;                                     \
      tests_failed++;                                                          \
    }                                                                          \
  } while (0)
