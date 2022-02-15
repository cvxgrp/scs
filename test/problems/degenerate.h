#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

/* test degenerate cones */
static const char *degenerate(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag, success;
  scs_float tpobj, tdobj, perr, derr;
  const char *fail;

  /* data */
  scs_float Ax[] = {-10., -1., 1., -1.};
  scs_int Ai[] = {1, 2, 1, 3};
  scs_int Ap[] = {0, 2, 4};

  scs_float Px[] = {0.02, 2.};
  scs_int Pi[] = {0, 1};
  scs_int Pp[] = {0, 1, 2};

  scs_float b[] = {1., -1., 2., -0.5};
  scs_float c[] = {1., 2.};

  scs_int m = 4;
  scs_int n = 2;

  /* used later: */
  scs_int sq[] = {1};

  /* end data */

  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));

  d->A->m = m;
  d->A->n = n;

  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  d->P->m = n;
  d->P->n = n;

  d->P->x = Px;
  d->P->i = Pi;
  d->P->p = Pp;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  /* positive orthants */
  k->l = 4;
  exitflag = scs(d, k, stgs, sol, &info);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) {
    return fail;
  }
  mu_assert("bsize: SCS failed to produce outputflag SCS_SOLVED",
            exitflag == SCS_SOLVED);

  tpobj = info.pobj;
  tdobj = info.dobj;

  /* degenerate box cone */
  k->bsize = 1;
  k->l = 3;
  exitflag = scs(d, k, stgs, sol, &info);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) {
    return fail;
  }
  perr = info.pobj - tpobj;
  derr = info.dobj - tdobj;
  success = ABS(perr) < 1e-8 && ABS(derr) < 1e-8 && exitflag == SCS_SOLVED;
  mu_assert("degenerate box cone failure", success);

  /* degenerate SOC cone */
  k->bsize = 0;
  k->q = sq;
  k->qsize = 1;
  k->l = 3;
  exitflag = scs(d, k, stgs, sol, &info);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) {
    return fail;
  }
  perr = info.pobj - tpobj;
  derr = info.dobj - tdobj;
  success = ABS(perr) < 1e-8 && ABS(derr) < 1e-8 && exitflag == SCS_SOLVED;
  mu_assert("degenerate SOC cone failure", success);

  /* degenerate PSD cone */
  k->q = 0;
  k->qsize = 0;
  k->s = sq;
  k->ssize = 1;
  k->l = 3;
  exitflag = scs(d, k, stgs, sol, &info);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) {
    return fail;
  }
  perr = info.pobj - tpobj;
  derr = info.dobj - tdobj;
  success = ABS(perr) < 1e-8 && ABS(derr) < 1e-8 && exitflag == SCS_SOLVED;
  mu_assert("degenerate PSD cone failure", success);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
