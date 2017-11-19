#include "scs.h"
#include "constants.h"
#include "../minunit.h"
#include "../problem_utils.h"

static const char *small_lp() {
  ScsCone *k = scs_calloc(1, sizeof(ScsCone));
  ScsData *d = scs_calloc(1, sizeof(ScsData));
  ScsSolution *sol = scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_float p_f = 0.1;
  int seed = 1234;
  scs_int n = 100;
  scs_int m = 300;
  scs_int col_nnz = (scs_int)ceil(sqrt(n));
  scs_int nnz = n * col_nnz;
  d->stgs = scs_calloc(1, sizeof(ScsSettings));
  k->f = (scs_int)floor(m * p_f);
  k->l = m - k->f;
  scs_int exitflag;

  srand(seed);

  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol);
  set_default_scs_settings(d);

  exitflag = scs(d, k, sol, &info);

  free_data(d, k);
  free_sol(sol);
  free_sol(opt_sol);
  mu_assert("small_lp: SCS failed to produce outputflag SCS_SOLVED",
            exitflag == SCS_SOLVED);
  return 0;
}
