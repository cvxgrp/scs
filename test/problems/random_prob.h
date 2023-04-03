#include "glbopts.h"
#include "scs.h"
#include "problems/prob_from_data_file.h"

static const char *random_prob(void) {
  scs_float OPT = 5.751458006385587;
  return _test_prob_from_data("test/problems/random_prob", OPT);
}
