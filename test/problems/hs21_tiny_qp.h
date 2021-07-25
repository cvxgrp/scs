#include "matrix.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "scs.h"
#include "util.h"


static const char *hs21_tiny_qp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char * fail;

  /* data */
  scs_float Ax[] = {-10.,  -1.,   1.,  -1.};
  scs_int Ai[] = {1, 2, 1, 3};
  scs_int Ap[] = {0, 2, 4};

  scs_float Px[] = {0.02, 2.};
  scs_int Pi[] = {0, 1};
  scs_int Pp[] = {0, 1, 2};

  scs_float b[] = {1.0, 0.0, 0.0};
  scs_float c[] ={0., 0.};

  scs_int m = 4;
  scs_int n = 2;

  scs_float bl[] = {10.0, 2.0, -50.0};
  scs_float bu[] = {1e+20, 50.0, 50.0};
  scs_int bsize = 3;

  scs_float opt = 0.04000000000000625;
  /* end data */

  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;

  d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
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

  k->bsize = bsize;
  k->bl = bl;
  k->bu = bu;

  SCS(set_default_settings)(d->stgs);
  d->stgs->eps_abs = 1e-6;
  d->stgs->eps_rel = 1e-6;
  d->stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("hs21_tiny_qp: SCS failed to produce outputflag SCS_SOLVED",
            success);
  fail = verify_solution_correct(d, k, &info, sol, exitflag);

  /* test warm-starting */
  d->stgs->warm_start = 1;
  exitflag = scs(d, k, sol, &info);
  /* 25 iters should be enough if warm-started */
  mu_assert("hs21_tiny_qp: warm-start failure", info.iter <= 25);
  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("hs21_tiny_qp: SCS failed to produce outputflag SCS_SOLVED",
            success);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(d->stgs);
  scs_free(d);
  return fail;
}
