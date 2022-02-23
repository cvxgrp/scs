#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "rw.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *rob_gauss_cov_est(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success, read_status;
  const char *fail;

  /* data */
  scs_float Ax[] = {-1.0,
                    -0.25035003149,
                    -0.214684983694,
                    -0.0013696512491,
                    -1.0,
                    -1.0,
                    0.0941450010301,
                    -0.0929840802212,
                    0.0143905024814,
                    -1.0,
                    -1.41421356237,
                    -1.0,
                    0.0941450010301,
                    -0.0929840802212,
                    0.0143905024814,
                    1.0,
                    -1.0,
                    -0.0354035554388,
                    -0.0402731435884,
                    -0.151196563215,
                    -1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    1.0,
                    -1.0,
                    1.0,
                    1.0,
                    -1.41421356237,
                    1.0,
                    -1.41421356237,
                    1.0,
                    -1.41421356237,
                    1.0,
                    -1.0,
                    1.0,
                    -1.41421356237,
                    1.0,
                    -1.41421356237,
                    1.0,
                    -1.0,
                    1.0,
                    1.0,
                    -1.41421356237,
                    1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    1.0,
                    -1.0,
                    -1.0,
                    1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0,
                    -1.0};
  scs_int Ai[] = {11, 15, 16, 17, 36, 12, 15, 16, 17, 18, 37, 13, 15, 16,
                  17, 18, 14, 15, 16, 17, 38, 15, 19, 16, 20, 17, 21, 19,
                  20, 21, 25, 19, 20, 21, 19, 22, 20, 23, 21, 24, 3,  26,
                  4,  5,  27, 7,  28, 9,  29, 6,  30, 8,  31, 10, 32, 11,
                  33, 12, 13, 34, 14, 35, 1,  7,  0,  8,  9,  2,  10, 1,
                  3,  41, 2,  6,  44, 25, 39, 25, 42};
  scs_int Ap[] = {0,  5,  11, 16, 21, 23, 25, 27, 31, 34, 36,
                  38, 40, 42, 45, 47, 49, 51, 53, 55, 57, 60,
                  62, 64, 66, 67, 69, 72, 75, 77, 79};
  scs_float b[] = {-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
                   -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
                   -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
                   -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
                   -0.0, -0.0, -0.0, -0.0, 1.0,  -0.0, -0.0, 1.0,  -0.0};
  scs_float c[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0,
                   1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  scs_int m = 45;
  scs_int n = 30;

  scs_int z = 19;
  scs_int l = 7;
  scs_int *q = SCS_NULL;
  scs_int qsize = 0;
  scs_int s[] = {4, 2};
  scs_int ssize = 2;
  scs_int ep = 2;
  scs_int ed = 0;
  scs_float *p = SCS_NULL;
  scs_int psize = 0;

  scs_float opt = -4.8912;
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

  scs_set_default_settings(stgs);
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;
  stgs->eps_infeas = 1e-9;

  exitflag = scs(d, k, stgs, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - opt;
  derr = -SCS(dot)(d->b, sol->y, d->m) - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("rob_gauss_cov_est: SCS failed to produce outputflag SCS_SOLVED",
            success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
  if (fail)
    return fail;

  /* test warm-starting */
  stgs->warm_start = 1;
  exitflag = scs(d, k, stgs, sol, &info);
  /* 100 iters should be enough if warm-started */
  mu_assert("rob_gauss_cov_est: warm-start failure", info.iter <= 100);
  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("rob_gauss_cov_est: SCS failed to produce outputflag SCS_SOLVED",
            success);

  /* test read / write */
  stgs->write_data_filename = "rob_gauss_cov_est";
  stgs->max_iters = 1;

  exitflag = scs(d, k, stgs, sol, &info);

  /* kill data */
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);

  read_status = SCS(read_data)("rob_gauss_cov_est", &d, &k, &stgs);

  if (read_status < 0) {
    return "Data read failure, exit.\n";
  }

  stgs->max_iters = 1000;
  /* solve with read data */
  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("rob_gauss_cov_est_rw: SCS failed to produce outputflag SCS_SOLVED",
            success);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  /* test other settings */
  stgs->warm_start = 0;
  stgs->normalize = 0;
  stgs->adaptive_scale = 0;
  stgs->acceleration_lookback = 10;
  stgs->acceleration_interval = 10;
  stgs->write_data_filename = SCS_NULL;
  stgs->max_iters = 1000;
  stgs->log_csv_filename = "rob_gauss_cov_est.csv";

  exitflag = scs(d, k, stgs, sol, &info);

  perr = info.pobj - opt;
  derr = info.dobj - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  mu_assert("rob_gauss_cov_est_rw: SCS failed to produce outputflag SCS_SOLVED",
            success);

  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);

  return fail;
}
