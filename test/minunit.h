/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */

/* Simple Macros for testing */

/* Failures are recorded rather than aborting the run, so that a single run
 * reports every failing test. The runner defines these. */
#define MU_MAX_FAILED_TESTS 256
extern int tests_run;
extern int tests_failed;
extern int tests_skipped;
extern const char *failed_tests[MU_MAX_FAILED_TESTS];
extern const char *skipped_tests[MU_MAX_FAILED_TESTS];
/* Sentinel returned by a test that is compiled out of the current build (see
 * _SKIP in run_tests.c). Compared by address, so it can never collide with a
 * real failure message. */
extern const char *const mu_skipped;
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
    if (message == mu_skipped) {                                               \
      scs_printf("skipped\n");                                                 \
      if (tests_skipped < MU_MAX_FAILED_TESTS)                                 \
        skipped_tests[tests_skipped] = name;                                   \
      tests_skipped++;                                                         \
    } else {                                                                   \
      tests_run++;                                                             \
      if (message) {                                                           \
        scs_printf("FAILED: %s: %s\n", name, message);                         \
        if (tests_failed < MU_MAX_FAILED_TESTS)                                \
          failed_tests[tests_failed] = name;                                   \
        tests_failed++;                                                        \
      }                                                                        \
    }                                                                          \
  } while (0)
