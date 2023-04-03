#include "glbopts.h"
#include "scs.h"
#include "problems/prob_from_data_file.h"

static const char *max_ent(void) {
  scs_float OPT = -6.0670878; /* from SCS so approximate */
  return _test_prob_from_data("test/problems/max_ent", OPT);
}