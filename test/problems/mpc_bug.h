#include "glbopts.h"
#include "problems/test_prob_from_data_file.h"
#include "scs.h"

static const char *mpc_bug(void) {
  const char *fail;
  scs_float OPT1 = -0.473957794500; /* from scs */
  scs_float OPT2 = -0.029336830816; /* from scs */
  scs_float OPT3 = -0.002215217478; /* from scs */
  fail = _test_prob_from_data("test/problems/mpc_bug1", OPT1);
  if (fail) {
    return fail;
  }
  fail = _test_prob_from_data("test/problems/mpc_bug2", OPT2);
  if (fail) {
    return fail;
  }
  return _test_prob_from_data("test/problems/mpc_bug3", OPT3);
}
