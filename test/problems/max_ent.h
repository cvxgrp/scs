#include "glbopts.h"
<<<<<<< HEAD
#include "problems/prob_from_data_file.h"
#include "scs.h"

static const char *max_ent(void) {
  scs_float OPT = -6.0670878; /* from SCS so approximate */
=======
#include "scs.h"
#include "problems/test_prob_from_data_file.h"

static const char *max_ent(void) {
  scs_float OPT = -6.067087663361563; /* from ecos */
>>>>>>> master
  return _test_prob_from_data("test/problems/max_ent", OPT);
}
