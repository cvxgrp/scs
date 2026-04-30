#include "glbopts.h"
#include "cones.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/* Shared minimal problem setup for validation tests. */
static void _setup_validation_prob(ScsData *d, ScsCone *k, ScsSolution *opt,
                                   int seed) {
  scs_float p_f = 0.1;
  scs_int n = 1, m = 3;
  scs_int col_nnz = 1, nnz = 1;
  k->z = (scs_int)(m * p_f);
  k->l = m - k->z;
  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt, seed);
}

static const char *test_validation(void) {
  scs_printf("Testing that SCS handles bad inputs correctly:\n");

  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol, *opt_sol;
  ScsInfo info = {0};
  scs_int exitflag;

#define VALIDATION_SETUP()                                               \
  do {                                                                   \
    k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));                       \
    d = (ScsData *)scs_calloc(1, sizeof(ScsData));                       \
    stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));            \
    sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));             \
    opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));         \
    _setup_validation_prob(d, k, opt_sol, 1234);                         \
    scs_set_default_settings(stgs);                                      \
  } while (0)

#define VALIDATION_CLEANUP()                                             \
  do {                                                                   \
    SCS(free_data)(d);                                                   \
    SCS(free_cone)(k);                                                   \
    SCS(free_sol)(sol);                                                  \
    SCS(free_sol)(opt_sol);                                              \
    scs_free(stgs);                                                      \
  } while (0)

  /* eps_abs < 0 */
  VALIDATION_SETUP();
  stgs->eps_abs = -1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: eps_abs < 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* eps_abs NaN */
  VALIDATION_SETUP();
  stgs->eps_abs = NAN;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: eps_abs NaN should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* eps_rel < 0 */
  VALIDATION_SETUP();
  stgs->eps_rel = -1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: eps_rel < 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* eps_infeas < 0 */
  VALIDATION_SETUP();
  stgs->eps_infeas = -1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: eps_infeas < 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* max_iters <= 0 */
  VALIDATION_SETUP();
  stgs->max_iters = 0;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: max_iters <= 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* alpha <= 0 */
  VALIDATION_SETUP();
  stgs->alpha = 0.0;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: alpha <= 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* alpha >= 2 */
  VALIDATION_SETUP();
  stgs->alpha = 2.0;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: alpha >= 2 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* rho_x <= 0 */
  VALIDATION_SETUP();
  stgs->rho_x = 0.0;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: rho_x <= 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* scale <= 0 */
  VALIDATION_SETUP();
  stgs->scale = 0.0;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: scale <= 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* time_limit_secs NaN */
  VALIDATION_SETUP();
  stgs->time_limit_secs = NAN;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: time_limit_secs NaN should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* acceleration_interval <= 0 */
  VALIDATION_SETUP();
  stgs->acceleration_interval = 0;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: acceleration_interval <= 0 should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* --- Cone validation tests --- */

  /* cone dimension mismatch: total cone dims != m */
  VALIDATION_SETUP();
  k->l += 99; /* inflate l so total dims > m */
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: cone dims != m should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* negative zero-cone dimension */
  VALIDATION_SETUP();
  k->l = 0;
  k->z = -1;
  d->m = -1; /* match total so cone-dims == m check passes */
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: z < 0 should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* box cone bl > bu */
  {
    scs_float bl_bad[] = {5.0};
    scs_float bu_bad[] = {1.0};
    VALIDATION_SETUP();
    k->l = 0;
    k->z = 0;
    k->bsize = 2; /* 1 t-var + 1 bounded var */
    k->bl = bl_bad;
    k->bu = bu_bad;
    d->m = 2;
    exitflag = scs(d, k, stgs, sol, &info);
    mu_assert("validation: bl > bu should fail", exitflag == SCS_FAILED);
    k->bl = SCS_NULL; /* prevent double-free */
    k->bu = SCS_NULL;
    VALIDATION_CLEANUP();
  }

  /* box cone one-sided infinite bounds are valid */
  {
    ScsData d_box = {0};
    ScsCone k_box = {0};
    scs_float bl_inf[] = {-INFINITY, -1.0};
    scs_float bu_inf[] = {1.0, INFINITY};
    d_box.m = 3;
    k_box.bsize = 3; /* 1 t-var + 2 bounded vars */
    k_box.bl = bl_inf;
    k_box.bu = bu_inf;
    mu_assert("validation: box infinite bounds should pass",
              SCS(validate_cones)(&d_box, &k_box) == 0);
  }

  /* box cone NaN bound */
  {
    ScsData d_box = {0};
    ScsCone k_box = {0};
    scs_float bl_nan[] = {NAN};
    scs_float bu_nan[] = {1.0};
    d_box.m = 2;
    k_box.bsize = 2; /* 1 t-var + 1 bounded var */
    k_box.bl = bl_nan;
    k_box.bu = bu_nan;
    mu_assert("validation: box NaN bound should fail",
              SCS(validate_cones)(&d_box, &k_box) < 0);
  }

  /* box cone infinities must point in the valid bound direction */
  {
    ScsData d_box = {0};
    ScsCone k_box = {0};
    scs_float bl_bad_inf[] = {INFINITY};
    scs_float bu_bad_inf[] = {INFINITY};
    d_box.m = 2;
    k_box.bsize = 2; /* 1 t-var + 1 bounded var */
    k_box.bl = bl_bad_inf;
    k_box.bu = bu_bad_inf;
    mu_assert("validation: box +inf lower bound should fail",
              SCS(validate_cones)(&d_box, &k_box) < 0);
  }

  /* power cone p out of range */
  {
    scs_float p_bad[] = {2.0}; /* must be in [-1, 1] */
    VALIDATION_SETUP();
    k->l = 0;
    k->z = 0;
    k->psize = 1;
    k->p = p_bad;
    d->m = 3;
    exitflag = scs(d, k, stgs, sol, &info);
    mu_assert("validation: p > 1 should fail", exitflag == SCS_FAILED);
    k->p = SCS_NULL; /* prevent double-free */
    VALIDATION_CLEANUP();
  }

  /* soc cone array missing */
  VALIDATION_SETUP();
  k->l = 0;
  k->z = 0;
  k->qsize = 1;
  k->q = SCS_NULL;
  d->m = 1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: missing soc cone array should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* missing settings */
  VALIDATION_SETUP();
  exitflag = scs(d, k, SCS_NULL, sol, &info);
  mu_assert("validation: missing settings should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* missing A */
  VALIDATION_SETUP();
  SCS(free_scs_matrix)(d->A);
  d->A = SCS_NULL;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: missing A should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* A dimensions inconsistent with ScsData */
  VALIDATION_SETUP();
  d->A->m = d->m + 1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: A dims inconsistent with data should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* A has negative row index */
  VALIDATION_SETUP();
  d->A->i[0] = -1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: A negative row index should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* A column pointers must be nondecreasing */
  VALIDATION_SETUP();
  d->A->p[0] = 1;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: A invalid column pointer should fail",
            exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

  /* A has non-finite data */
  VALIDATION_SETUP();
  d->A->x[0] = NAN;
  exitflag = scs(d, k, stgs, sol, &info);
  mu_assert("validation: A NaN entry should fail", exitflag == SCS_FAILED);
  VALIDATION_CLEANUP();

#undef VALIDATION_SETUP
#undef VALIDATION_CLEANUP

  return 0;
}
