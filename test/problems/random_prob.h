#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "rw.h"
#include "scs.h"
#include "util.h"

static const char *random_prob(void) {
  scs_int read_status;
  ScsData *d;
  ScsCone *k;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  scs_float OPT = 5.751458006385587;

  read_status = SCS(read_data)("test/problems/random_prob", &d, &k, &stgs);

  if (read_status < 0) {
    return "Data read failure, exit.\n";
  }

  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  exitflag = scs(d, k, stgs, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - OPT;
  derr = -SCS(dot)(d->b, sol->y, d->m) - OPT;
  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("random_prob: SCS failed to produce SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);

  return fail;
}
