#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *darwinn2(void) {
  scs_printf("\nDarwinn2 Test Begin\n");
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  /* data */
  scs_float Ax[] = {
      -1,          -1, -1,          -1, -1,          1,  1,           -1,
      -1,          -1, -1,          -1, 1,           1,  -1,          -1,
      -1,          -1, -1,          -1, -1,          -1, -1,          -1,
      6.20355e-10, -1, 8.17436e-06, -1, 6.05508e-06, -1, 0.000341506, -1,
      0.000178019, -1, 6.05508e-07, -1, 9.08262e-06, -1, 6.66059e-06, -1,
      0.000343323, -1, 0.93196,     -1, 0.0729925,   -1, 1.26694,     -1,
      66.9777,     -1, -1};
  scs_int Ai[] = {5,  8,  11, 23, 26, 35, 44, 6,  8,  14, 23, 29, 41,
                  44, 7,  8,  11, 14, 17, 20, 23, 26, 29, 32, 0,  10,
                  0,  13, 0,  16, 0,  19, 1,  22, 2,  25, 2,  28, 2,
                  31, 2,  34, 3,  37, 3,  40, 3,  43, 3,  46, 3};
  scs_int Ap[] = {0,  7,  14, 24, 26, 28, 30, 32, 34,
                  36, 38, 40, 42, 44, 46, 48, 50, 51};

  scs_float b[] = {/* cone_l */
                   1, 1, 1, 0,
                   /* cone_box */
                   1, 0, 0, 0,
                   /* cone_ep */
                   0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
                   0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};
  scs_float c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  scs_int m = 47;
  scs_int n = 17;

  scs_int cone_l = 4;
  scs_float bl[] = {0.0, 0.0, 0.0};
  scs_float bu[] = {6.23832, 6.51026, 2.07944};
  scs_int cone_bsize = 4;
  scs_int cone_ep = 13;

  scs_float opt = 0.07689269067639115;
  /* end data */

  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;

  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));

  d->A->m = m;
  d->A->n = n;

  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  k->l = cone_l;
  k->bl = bl;
  k->bu = bu;
  k->bsize = cone_bsize;
  k->ep = cone_ep;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-7;
  stgs->eps_rel = 1e-7;
  stgs->eps_infeas = 1e-9;

  stgs->normalize = 0;
  stgs->max_iters = 1000000;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("darwinn2: SCS failed to produce outputflag SCS_SOLVED", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  scs_printf("\nDarwinn2 Test End\n");
  return fail;
}
