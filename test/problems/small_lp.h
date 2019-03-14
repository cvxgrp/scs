#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

static const char *small_lp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_float p_f = 0.1;
  int seed = 1234;
  scs_int n = 100;
  scs_int m = 300;
  scs_int col_nnz = (scs_int)ceil(sqrt(n));
  scs_int nnz = n * col_nnz;
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;

  d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  k->f = (scs_int)floor(m * p_f);
  k->l = m - k->f;

  srand(seed);

  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol);
  SCS(set_default_settings)(d);

  /* this test fails with the default choice of 10 */
  d->stgs->acceleration_lookback = 20;

  exitflag = scs(d, k, sol, &info);
  perr = SCS(dot)(d->c, sol->x, d->n) - SCS(dot)(d->c, opt_sol->x, d->n);
  derr = -SCS(dot)(d->b, sol->y, d->m) + SCS(dot)(d->b, opt_sol->y, d->m);
  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  SCS(free_data)(d, k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  mu_assert("small_lp: SCS failed to produce outputflag SCS_SOLVED", success);
  return 0;
}
