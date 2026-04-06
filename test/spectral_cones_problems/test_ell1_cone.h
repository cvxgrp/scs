#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

/*
 * Basis pursuit: minimize ||x||_1 subject to Ax = b
 *
 * Conic form:
 *   minimize t
 *   subject to:
 *     x1 + 2*x2     = 1   (equality)
 *          x2  + x3  = 1   (equality)
 *     (t, x1, x2, x3) in ell1(3)
 *
 * Optimal: x = (0, 0.5, 0.5), t = 1.0
 */
static const char *test_ell1_cone(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;

  /* A matrix (6 x 4, CSC):
   *  Row 0: [0,  1,  2,  0]   (zero)
   *  Row 1: [0,  0,  1,  1]   (zero)
   *  Row 2: [-1, 0,  0,  0]   (ell1 t)
   *  Row 3: [0, -1,  0,  0]   (ell1 x1)
   *  Row 4: [0,  0, -1,  0]   (ell1 x2)
   *  Row 5: [0,  0,  0, -1]   (ell1 x3)
   */
  scs_float Ax[] = {-1., 1., -1., 2., 1., -1., 1., -1.};
  scs_int Ai[] = {2, 0, 3, 0, 1, 4, 1, 5};
  scs_int Ap[] = {0, 1, 3, 6, 8};
  scs_float b[] = {1., 1., 0., 0., 0., 0.};
  scs_float c[] = {1., 0., 0., 0.};

  scs_int m = 6;
  scs_int n = 4;

  scs_int ell1[] = {3};
  scs_int ell1_size = 1;

  scs_int *d_array = SCS_NULL;
  scs_int dsize = 0;
  scs_int *nuc_m = SCS_NULL;
  scs_int *nuc_n = SCS_NULL;
  scs_int nucsize = 0;
  scs_int *sl_n = SCS_NULL;
  scs_int *sl_k = SCS_NULL;
  scs_int sl_size = 0;

  scs_float opt = 1.0;

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

  k->z = 2;
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
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - opt;
  derr = -SCS(dot)(d->b, sol->y, d->m) - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("test_ell1_cone: SCS failed to produce outputflag SCS_SOLVED",
            success);

  mu_assert("test_ell1_cone: primal feas error", ABS(info.res_pri) < 1e-6);
  mu_assert("test_ell1_cone: dual feas error", ABS(info.res_dual) < 1e-6);
  mu_assert("test_ell1_cone: duality gap error", ABS(info.gap) < 1e-6);

  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  SCS(free_sol)(sol);

  return 0;
}
