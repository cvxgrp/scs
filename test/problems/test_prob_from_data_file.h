#ifndef _SCS_FILE_TEST_CHASSIS
#define _SCS_FILE_TEST_CHASSIS

#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "rw.h"
#include "scs.h"
#include "util.h"

static const char *_test_prob_from_data(const char *file, scs_float OPT) {
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
  scs_float xt_p_x;
  scs_float *px = SCS_NULL;

  read_status = SCS(read_data)(file, &d, &k, &stgs);

  if (read_status < 0) {
    return "Data read failure, exit.\n";
  }

  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  /* Force verbosity for the test */
  stgs->verbose = 1;

  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  exitflag = scs(d, k, stgs, sol, &info);

  if (d->P) {
    /* px = Px */
    px = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
    memset(px, 0, d->n * sizeof(scs_float));
    SCS(accum_by_p)(d->P, sol->x, px);
    xt_p_x = SCS(dot)(px, sol->x, d->n);
  } else {
    xt_p_x = 0.;
  }

  perr = 0.5 * xt_p_x + SCS(dot)(d->c, sol->x, d->n) - OPT;
  derr = -0.5 * xt_p_x - SCS(dot)(d->b, sol->y, d->m) - OPT;
  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;
  if (!success) {
    scs_printf("%s: FAILED\n", file);
  }
  mu_assert(file, success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);
  if (px) {
    scs_free(px);
  }
  if (fail) {
    scs_printf("%s: FAILED\n", file);
  }
  return fail;
}

#endif
