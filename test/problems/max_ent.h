#include "glbopts.h"
#include "scs.h"
#include "problems/test_prob_from_data_file.h"

static const char *max_ent(void) {
  scs_float OPT = -6.067087663361563; /* from ecos */
  return _test_prob_from_data("test/problems/max_ent", OPT);
}
