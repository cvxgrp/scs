#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

static const char *test_validation(void) {
  scs_printf("Testing that SCS handles bad inputs correctly:\n");

  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_float p_f = 0.1;
  int seed = 1234;
  scs_int n = 1;
  scs_int m = 3;
  scs_int col_nnz = (scs_int)ceil(sqrt(n));
  scs_int nnz = n * col_nnz;
  scs_int exitflag;

  k->z = (scs_int)floor(m * p_f);
  k->l = m - k->z;

  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed);
  scs_set_default_settings(stgs);

  /* TODO test more failure modes */
  stgs->eps_abs = -1;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("test_fails: SCS failed to produce outputflag SCS_FAILED",
            exitflag == SCS_FAILED);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);

  return 0;
}
