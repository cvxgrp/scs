#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * cvxgrp/scs#499: a barely infeasible feasibility problem (c = 0, no P), one
 * of the subproblems of a DQCP bisection (nested_call011 in the issue's
 * reproducer, same data and the settings CVXR passes). m = 12, n = 4, three
 * nonnegative rows and three second-order cones of size 3.
 *
 * With c = 0 and P = 0 the dual residual is A'y and so is the only term it
 * is measured against: its relative size is 1 on every iteration. A scale
 * update that balances relative residuals therefore always reads the dual
 * as behind and walks the scale down to its floor, where progress toward
 * the certificate stalls. 3.3.1 needed 84,475 iterations here (3.2.7:
 * 5,500); with the residuals themselves balanced whenever a side has no
 * correct digits it takes about 1,250.
 *
 * The counts move with the platform and the BLAS, so the test only asks for
 * the certificate within 20,000 iterations.
 */
static const char *issue_499(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  scs_int Ap[] = {0, 5, 8, 10, 14};
  scs_int Ai[] = {0, 2, 3, 4, 8, 1, 3, 4, 1, 11, 6, 7, 9, 10};
  scs_float Ax[] = {-1., -1.,   -1., -1., -2., 0.05, -1.,
                    1.,  -0.35, -2., -1., -1., 1.,   1.};
  scs_float b[] = {0., 0., -0.9934259201488199, 0., 0., 2., 1., -1., 0., 2.,
                   0., 0.};
  scs_float c[] = {0., 0., 0., 0.};
  scs_int q[] = {3, 3, 3};

  scs_int m = 12, n = 4;

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

  k->l = 3;
  k->q = q;
  k->qsize = 3;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-5;
  stgs->eps_rel = 1e-5;
  stgs->acceleration_lookback = 10;
  stgs->acceleration_interval = 10;
  stgs->max_iters = 20000;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("issue_499: expected SCS_INFEASIBLE within 20,000 iterations",
            exitflag == SCS_INFEASIBLE);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
