#include "glbopts.h"
#include "problems/test_prob_from_data_file.h"
#include "scs.h"

static const char *mpc_bug(void) {
  scs_float OPT = -0.473957794500; /* from scs */
  return _test_prob_from_data("test/problems/mpc_bug", OPT);
}
