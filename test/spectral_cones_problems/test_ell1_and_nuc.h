#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

/*
 * Combined ell1 + nuclear norm test.
 * Two independent blocks sharing one SCS call to exercise cone ordering.
 *
 * Block 1 (ell1):
 *   minimize t_ell s.t. a1 + a2 = 1, (t_ell, a1, a2) in ell1(2)
 *   Optimal: t_ell = 1
 *
 * Block 2 (nuclear norm):
 *   minimize t_nuc s.t. M = [[0.6, 0], [0.8, 0]], (t_nuc, M) in nuc(2,2)
 *   M has singular values {1, 0}, so optimal t_nuc = 1
 *
 * Variables (n=8): t_ell, a1, a2, t_nuc, M_11, M_21, M_12, M_22
 * Objective: minimize t_ell + t_nuc, optimal = 2.0
 *
 * Cone ordering: z=5, nuc(2,2)=5 rows, ell1(2)=3 rows, total m=13
 */
static const char *test_ell1_and_nuc(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;

  /* A matrix (13 x 8, CSC):
   *  Row  0: [0, 1, 1, 0, 0, 0, 0, 0]   z: a1 + a2 = 1
   *  Row  1: [0, 0, 0, 0, 1, 0, 0, 0]   z: M_11 = 0.6
   *  Row  2: [0, 0, 0, 0, 0, 1, 0, 0]   z: M_21 = 0.8
   *  Row  3: [0, 0, 0, 0, 0, 0, 1, 0]   z: M_12 = 0
   *  Row  4: [0, 0, 0, 0, 0, 0, 0, 1]   z: M_22 = 0
   *  Row  5: [0, 0, 0,-1, 0, 0, 0, 0]   nuc: t_nuc
   *  Row  6: [0, 0, 0, 0,-1, 0, 0, 0]   nuc: M_11
   *  Row  7: [0, 0, 0, 0, 0,-1, 0, 0]   nuc: M_21
   *  Row  8: [0, 0, 0, 0, 0, 0,-1, 0]   nuc: M_12
   *  Row  9: [0, 0, 0, 0, 0, 0, 0,-1]   nuc: M_22
   *  Row 10: [-1,0, 0, 0, 0, 0, 0, 0]   ell1: t_ell
   *  Row 11: [0,-1, 0, 0, 0, 0, 0, 0]   ell1: a1
   *  Row 12: [0, 0,-1, 0, 0, 0, 0, 0]   ell1: a2
   */
  scs_float Ax[] = {-1., 1., -1., 1., -1., -1., 1., -1., 1., -1.,
                    1., -1., 1., -1.};
  scs_int Ai[] = {10, 0, 11, 0, 12, 5, 1, 6, 2, 7, 3, 8, 4, 9};
  scs_int Ap[] = {0, 1, 3, 5, 6, 8, 10, 12, 14};
  scs_float b[] = {1., 0.6, 0.8, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  scs_float c[] = {1., 0., 0., 1., 0., 0., 0., 0.};

  scs_int m = 13;
  scs_int n = 8;

  scs_int nuc_m[] = {2};
  scs_int nuc_n[] = {2};
  scs_int nucsize = 1;
  scs_int ell1[] = {2};
  scs_int ell1_size = 1;

  scs_int *d_array = SCS_NULL;
  scs_int dsize = 0;
  scs_int *sl_n = SCS_NULL;
  scs_int *sl_k = SCS_NULL;
  scs_int sl_size = 0;

  scs_float opt = 2.0;

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

  k->z = 5;
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

  mu_assert("test_ell1_and_nuc: SCS failed to produce outputflag SCS_SOLVED",
            success);

  mu_assert("test_ell1_and_nuc: primal feas error", ABS(info.res_pri) < 1e-6);
  mu_assert("test_ell1_and_nuc: dual feas error", ABS(info.res_dual) < 1e-6);
  mu_assert("test_ell1_and_nuc: duality gap error", ABS(info.gap) < 1e-6);

  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  SCS(free_sol)(sol);

  return 0;
}
