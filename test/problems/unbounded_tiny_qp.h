#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *unbounded_tiny_qp(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  const char *fail;

  /* data */
  scs_float Ax[] = {
      -0.04101197, 1.68314201,  0.5543439,   -0.96712378, 0.19370457,
      -1.15055711, 1.13315795,  -1.05910693, 1.24512977,  -0.70999048,
      -0.89976326, -0.41373294, -0.73848186, 0.08663554,  -0.1681749,
      -0.42711619, -0.90501247, -0.490446,   -0.67124734, 1.67772257,
      0.39924394,  0.16330292,  0.55609205,  -1.22088238, -0.25891675,
      -3.07159984, -1.84102417, 1.5621635,   -1.13771529, 0.56067264,
      -0.0854747,  0.31024722,  -0.07437118, -0.20711534, -0.35241366,
      -0.98965142, -1.91488894, 1.01591507,  0.45387459,  1.43709968,
      -0.0482982,  -0.32447,    -0.91433399, 0.49750765,  0.09150015,
      0.69164184,  0.51064936,  -1.35009809, -1.35403213, 1.51897627};
  scs_int Ai[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6,
                  7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3,
                  4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  scs_int Ap[] = {0, 10, 20, 30, 40, 50};

  scs_float b[] = {1.4672744,   0.36190605, -0.5548082,  -0.35039932,
                   -1.12765224, 0.51012137, -0.24927975, -1.45270362,
                   -1.94080389, -0.0189713};
  scs_float c[] = {0.17547329, -0.11983635, -0.11791039, 0.12099476,
                   0.61931906};

  scs_int m = 10;
  scs_int n = 5;

  scs_int l = m;

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

  k->l = l;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-7;

  exitflag = scs(d, k, stgs, sol, &info);

  mu_assert("unbounded_tiny_qp: SCS failed to produce outputflag SCS_UNBOUNDED",
            exitflag == SCS_UNBOUNDED);

  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  SCS(free_cone)(k);
  scs_free(stgs);
  scs_free(d->A);
  scs_free(d);
  return fail;
}
