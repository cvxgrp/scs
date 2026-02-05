#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *partial_warm_start(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_int cold_iters;
  scs_int i;
  scs_float perr, derr;
  scs_int success;
  const char *fail;

  /* data (same as hs21_tiny_qp) */
  scs_float Ax[] = {-10., -1., 1., -1.};
  scs_int Ai[] = {1, 2, 1, 3};
  scs_int Ap[] = {0, 2, 4};

  scs_float Px[] = {0.02, 2.};
  scs_int Pi[] = {0, 1};
  scs_int Pp[] = {0, 1, 2};

  scs_float b[] = {1., 0., 0., 0.};
  scs_float c[] = {0., 0.};

  scs_int m = 4;
  scs_int n = 2;

  scs_float bl[] = {10.0, 2.0, -50.0};
  scs_float bu[] = {1e+20, 50.0, 50.0};
  scs_int bsize = 4;

  scs_float opt = 0.04000000000000625;
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

  k->bsize = bsize;
  k->bl = bl;
  k->bu = bu;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;
  stgs->eps_infeas = 0.;
  stgs->acceleration_lookback = 0; /* disable acceleration for consistent iters */

  /* Step 1: Cold solve to get reference solution and iteration count */
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start: cold solve failed", success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail) {
    SCS(free_sol)(sol);
    scs_free(d->A);
    scs_free(d->P);
    scs_free(k);
    scs_free(stgs);
    scs_free(d);
    return fail;
  }

  cold_iters = info.iter;
  scs_printf("partial_warm_start: cold solve took %li iters\n",
             (long)cold_iters);

  /* Step 2: Partial warm start - keep x, set s and y to NaN */
  for (i = 0; i < m; ++i) {
    sol->s[i] = NAN;
    sol->y[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start: partial warm start (x only) failed to solve",
            success);

  scs_printf("partial_warm_start: partial warm start (x only) took %li iters\n",
             (long)info.iter);
  mu_assert(
      "partial_warm_start: partial warm start should take fewer iters than "
      "cold start",
      info.iter < cold_iters);

  /* Step 3: All-NaN warm start - should still converge */
  for (i = 0; i < n; ++i) {
    sol->x[i] = NAN;
  }
  for (i = 0; i < m; ++i) {
    sol->s[i] = NAN;
    sol->y[i] = NAN;
  }

  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  success = ABS(perr) < 1e-3 && ABS(derr) < 1e-3 && exitflag == SCS_SOLVED;
  mu_assert("partial_warm_start: all-NaN warm start failed to solve", success);

  scs_printf("partial_warm_start: all-NaN warm start took %li iters\n",
             (long)info.iter);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
