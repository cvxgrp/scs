#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *infeasible_tiny_qp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  /* data */
  scs_float Ax[] = {
      4.51689976e-01,  3.06592046e-03,  5.17304192e-01,  -3.03038477e+00,
      -1.40509892e+00, 7.94277342e-04,  2.39454841e+00,  -7.60957360e-01,
      -1.18946302e+00, 3.98797701e-01,  1.01386914e+00,  -2.53921734e+00,
      -3.21202445e-01, 9.25735134e-01,  -2.54046934e-01, 1.28211442e-01,
      -1.65155072e-01, -4.53308401e-01, -4.66709068e-01, -2.24298562e-01,
      -4.92029627e-01, 8.05750411e-01,  -1.72920210e+00, -1.45633836e-01,
      6.39086786e-01,  1.20509649e-01,  4.19672104e-01,  -2.29274817e-01,
      3.30125838e-01,  4.12874191e-01,  1.05357823e+00,  2.10587878e+00,
      5.54934230e-01,  2.42617608e+00,  -8.33596918e-01, -6.83444334e-01,
      8.23856780e-01,  -4.14240741e-01, -7.24051659e-01, 1.06144329e-01,
      6.92857027e-01,  1.64822980e+00,  1.94603528e-01,  -7.58318366e-01,
      1.42948833e+00,  -2.49039902e-01, 2.15319490e-01,  -1.52651434e+00,
      -3.94761305e-01, -5.24305949e-01};
  scs_int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6,
                  7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3,
                  4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  scs_int Ap[] = {0, 10, 20, 30, 40, 50};

  scs_float Px[] = {0.40113649, 0.15897142, 0.37369516,  -0.05899464,
                    0.0772996,  0.37333677, -0.02994108, 0.04870077,
                    0.32221413, 0.28665235, -0.27383424, 0.04842465,
                    0.07832219, 0.0660133,  0.49993284};
  scs_int Pi[] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4};
  scs_int Pp[] = {0, 1, 3, 6, 10, 15};

  scs_float b[] = {-0.057279,   0.10451534, -0.10459346, -0.18452294,
                   0.18609658,  0.15795524, -0.00562746, -0.85852899,
                   -0.48271335, 0.61951655};
  scs_float c[] = {-0.16566256, -1.33116808, -0.1767858, -1.0940148,
                   1.15348983};

  scs_int m = 10;
  scs_int n = 5;

  scs_int l = m;

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

  k->l = l;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert(
      "infeasible_tiny_qp: SCS failed to produce outputflag SCS_INFEASIBLE",
      exitflag == SCS_INFEASIBLE);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
