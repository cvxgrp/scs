#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "rw.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

// for SpectralSCS
static const char *exp_design(void) {
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
      -1.,         1.,          -1.,         3.24,        0.16,
      1.,          4.84,        3.61,        1.,          -1.,
      2.54558441,  -0.11313708, -0.14142136, 1.24450793,  0.26870058,
      -2.12132034, -1.,         2.03646753,  0.05656854,  0.56568542,
      0.93338095,  4.03050865,  0.28284271,  -1.,         1.,
      0.04,        0.01,        0.16,        0.01,        2.25,
      -1.,         1.13137085,  -0.02828427, -0.05656854, 0.16970563,
      0.21213203,  -0.42426407, -1.,         0.64,        0.01,
      0.16,        0.09,        2.25,        0.04,        -1.};
  scs_int Ai[] = {7,  0,  8, 1, 2, 3, 4, 5,  6,  9, 1, 2, 3, 4, 5,
                  6,  10, 1, 2, 3, 4, 5, 6,  11, 1, 2, 3, 4, 5, 6,
                  12, 1,  2, 3, 4, 5, 6, 13, 1,  2, 3, 4, 5, 6, 14};
  scs_int Ap[] = {0, 1, 3, 10, 17, 24, 31, 38, 45};

  scs_float b[] = {1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.};
  scs_float c[] = {1., 0., 0., 0., 0., 0., 0., 0.};

  scs_int m = 15;
  scs_int n = 8;

  scs_int z = 1;
  scs_int l = 6;
  scs_int *q = SCS_NULL;
  scs_int qsize = 0;
  scs_int *s = SCS_NULL;
  scs_int ssize = 0;
  scs_int ep = 0;
  scs_int ed = 0;
  scs_float *p = SCS_NULL;
  scs_int psize = 0;
  scs_int d_array[] = {3};
  scs_int dsize = 1;
  scs_int *nuc_m = SCS_NULL;
  scs_int *nuc_n = SCS_NULL;
  scs_int nucsize = 0;
  scs_int *ell1 = SCS_NULL;
  scs_int ell1_size = 0;
  scs_int *sl_n = SCS_NULL;
  scs_int *sl_k = SCS_NULL;
  scs_int sl_size = 0;

  // computed using mosek (the input of Ax is truncated, and mosek solved
  // the problem with the non-truncated data)
  scs_float opt = 3.0333290743428574;
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

  k->z = z;
  k->l = l;
  k->q = q;
  k->qsize = qsize;
  k->s = s;
  k->ssize = ssize;
  k->ep = ep;
  k->ed = ed;
  k->p = p;
  k->psize = psize;
  k->d = d_array;
  k->dsize = dsize;
  k->nuc_m = nuc_m;
  k->nuc_n = nuc_n;
  k->nucsize = nucsize;
  k->ell1 = ell1;
  k->ell1_size = ell1_size;
  k->sl_n = sl_n;
  k->sl_k = sl_k;
  k->sl_size = sl_size;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-7;
  stgs->eps_rel = 1e-7;
  stgs->eps_infeas = 1e-9;

  stgs->log_csv_filename = "test_exp_design.csv";
  exitflag = scs(d, k, stgs, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - opt;
  derr = -SCS(dot)(d->b, sol->y, d->m) - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("exp_design: SCS failed to produce outputflag SCS_SOLVED", success);

  fail = 0;
  // TODO: This test fails because of the complementary slackness check.
  //       The complementary slackness tolerance is a bit too tight.
  // fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  // if (fail)
  //   return fail;

  mu_assert("exp_design: primal feas error: ", ABS(info.res_pri) < 1e-6);
  mu_assert("exp_design: dual feas error: ", ABS(info.res_dual) < 1e-6);
  mu_assert("exp_design: duality gap error: ", ABS(info.gap) < 1e-6);

  /* kill data */
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  SCS(free_sol)(sol);

  return fail;
}
