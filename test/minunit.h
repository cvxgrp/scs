/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */

/* Simple Macros for testing */
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
    if (message)                                                               \
      return message;                                                          \
  } while (0)
