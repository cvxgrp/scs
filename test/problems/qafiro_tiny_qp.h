#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *qafiro_tiny_qp(void) {
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
      -1.,    -1.06, -1.,    -0.301, -1.,    1.,     1.,     -1.,   1.,
      1.,     -1.,   1.,     -1.,    -1.,    -1.,    -1.06,  -1.,   -0.301,
      -1.,    -1.,   -1.06,  -1.,    -0.313, -1.,    -1.,    -0.96, -1.,
      -0.313, -1.,   -1.,    -0.86,  -1.,    -0.326, -1.,    1.,    -2.364,
      -1.,    1.,    -2.386, -1.,    1.,     -2.408, -1.,    1.,    -2.429,
      -1.,    1.,    -1.4,   -1.,    1.,     1.,     -1.,    1.,    -1.,
      -1.,    -0.43, -1.,    -1.,    -0.109, -1.,    1.,     1.,    -1.,
      1.,     1.,    -1.,    1.,     1.,     -1.,    1.,     -1.,   -1.,
      1.,     -0.43, -1.,    -0.109, -1.,    1.,     -0.43,  -1.,   -0.108,
      -1.,    1.,    -0.39,  -1.,    -0.108, -1.,    1.,     -0.37, -1.,
      -0.107, -1.,   1.,     -2.191, -1.,    1.,     -2.219, -1.,   1.,
      -2.249, -1.,   1.,     -2.279, -1.,    -1.,    -1.4,   -1.,   1.,
      1.,     -1.,   1.,     -1.,    -1.,    1.,     -1.};
  scs_int Ai[] = {0,  1,  16, 24, 28, 0,  15, 29, 0,  22, 30, 1,  26, 31, 4,
                  5,  12, 25, 32, 4,  5,  11, 25, 33, 4,  5,  9,  25, 34, 4,
                  5,  10, 25, 35, 12, 21, 36, 11, 21, 37, 9,  21, 38, 10, 21,
                  39, 4,  15, 40, 4,  23, 41, 5,  27, 42, 6,  7,  13, 22, 43,
                  7,  14, 44, 7,  24, 45, 7,  21, 46, 6,  26, 47, 2,  3,  17,
                  23, 48, 2,  3,  18, 23, 49, 2,  3,  19, 23, 50, 2,  3,  20,
                  23, 51, 17, 21, 52, 18, 21, 53, 19, 21, 54, 20, 21, 55, 2,
                  14, 56, 2,  25, 57, 3,  27, 58, 2,  59};
  scs_int Ap[] = {0,  5,  8,  11, 14, 19,  24,  29,  34,  37,  40,
                  43, 46, 49, 52, 55, 60,  63,  66,  69,  72,  77,
                  82, 87, 92, 95, 98, 101, 104, 107, 110, 113, 115};

  scs_float Px[] = {10., 1., 10., 1., 1., 10.};
  scs_int Pi[] = {0, 0, 1, 0, 1, 2};
  scs_int Pp[] = {0, 1, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};

  scs_float b[] = {
      0.00000000e+00, 0.00000000e+00, 4.40000000e+01, -2.22044605e-16,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00};
  scs_float c[] = {0., -0.4,  0., 0., 0., 0.,   0.,    0., 0., 0., 0.,
                   0., -0.32, 0., 0., 0., -0.6, 0.,    0., 0., 0., 0.,
                   0., 0.,    0., 0., 0., 0.,   -0.48, 0., 0., 10.};

  scs_int m = 60;
  scs_int n = 32;

  scs_float bl[] = {
      -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20,
      -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20, -1e+20,
      -1e+20, 0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      0.0,    0.0,    0.0,    0.0,    0.0,    0.0};
  scs_float bu[] = {0.0,   0.0,   0.0,   80.0,  500.0, 0.0,   0.0,   80.0,
                    500.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                    0.0,   310.0, 300.0, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20,
                    1e+20, 1e+20, 1e+20};
  scs_int bsize = 52;
  scs_int z = 8;

  scs_float opt = -1.5907818;
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
  k->z = z;

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-7;
  stgs->eps_rel = 1e-7;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("qafiro_tiny_qp: SCS failed to produce outputflag SCS_SOLVED",
            success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  /* test warm-starting */
  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);
  /* 25 iters should be enough if warm-started */
  mu_assert("qafiro_tiny_qp: warm-start failure", info.iter <= 100);
  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("qafiro_tiny_qp: SCS failed to produce outputflag SCS_SOLVED",
            success);

  /* test other settings */
  stgs->warm_start = 0;
  stgs->normalize = 0;
  stgs->adaptive_scale = 0;
  stgs->acceleration_lookback = 10;
  stgs->acceleration_interval = 10;
  stgs->log_csv_filename = "qafiro_tiny_qp.csv";

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("qafiro_tiny_qp: SCS failed to produce outputflag SCS_SOLVED",
            success);

  /* test other settings */
  stgs->warm_start = 0;
  stgs->normalize = 0;
  stgs->adaptive_scale = 0;
  stgs->acceleration_lookback = -10;
  stgs->acceleration_interval = 10;
  stgs->log_csv_filename = "qafiro_tiny_qp.csv";

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("qafiro_tiny_qp: SCS failed to produce outputflag SCS_SOLVED",
            success);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
