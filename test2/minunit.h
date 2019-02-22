/* Taken from http://www.jera.com/techinfo/jtns/jtn002.html */

/* Simple Macros for testing */
#define mu_assert(message, test) \
  do {                           \
    if (!(test)) return message; \
  } while (0)
#define mu_run_test(test)         \
  do {                            \
    const char *message = test(); \
    tests_run++;                  \
    if (message) return message;  \
  } while (0)
