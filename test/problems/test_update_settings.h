#include <stdio.h>
#include <string.h>

#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_work.h"
#include "util.h"

/*
 * Tests of scs_update_settings: changing the settings that only the solve
 * reads on an existing workspace, and refusing a change to a setting that
 * was used to set the workspace up.
 *
 * All of them use the random LP of small_lp (n = 100, m = 300), which takes
 * well over one iteration from a cold start, so a time limit of 1e-10 s
 * always stops it at the first convergence check (iteration 0), as in
 * test_time_limit_secs.
 */

#define UPD_N 100
#define UPD_M 300
#define UPD_SEED 1234

/* Builds the LP; the caller frees d, k and opt_sol. */
static void upd_gen_lp(ScsData **d, ScsCone **k, ScsSolution **opt_sol) {
  scs_int col_nnz = (scs_int)ceil(sqrt(UPD_N));
  *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  *opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  (*k)->z = (scs_int)floor(UPD_M * 0.1);
  (*k)->l = UPD_M - (*k)->z;
  (*d)->m = UPD_M;
  (*d)->n = UPD_N;
  gen_random_prob_data(UPD_N * col_nnz, col_nnz, *d, *k, *opt_sol, UPD_SEED);
}

static void upd_free(ScsData *d, ScsCone *k, ScsSolution *opt_sol,
                     ScsSolution *sol) {
  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(opt_sol);
  SCS(free_sol)(sol);
}

static int upd_hit_time_limit(scs_int exitflag, const ScsInfo *info) {
  return (exitflag == SCS_SOLVED_INACCURATE ||
          exitflag == SCS_INFEASIBLE_INACCURATE ||
          exitflag == SCS_UNBOUNDED_INACCURATE) &&
         strstr(info->status, "time_limit_secs") != SCS_NULL;
}

/*
 * The use case that motivated the call: set the time limit after setup.
 * A time limit set after scs_init takes effect in the next solve, stays in
 * force for the solves after it, and can be removed again; and one given
 * to scs_init can be removed before the first solve. The setup time that
 * ScsInfo reports is not changed by the call.
 */
static const char *test_update_settings_time_limit(void) {
  ScsData *d;
  ScsCone *k;
  ScsSolution *opt_sol;
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSettings stgs;
  ScsInfo info = {0};
  ScsWork *w;
  scs_int exitflag;
  scs_float setup_time;
  const char *fail = 0;

  upd_gen_lp(&d, &k, &opt_sol);
  scs_set_default_settings(&stgs);
  stgs.verbose = 1;

  w = scs_init(d, k, &stgs);
  mu_assert("update_settings_time_limit: scs_init failed", w);

  exitflag = scs_solve(w, sol, &info, 0);
  setup_time = info.setup_time;
  if (exitflag != SCS_SOLVED) {
    fail = "update_settings_time_limit: first solve (no limit) not solved";
    goto out;
  }

  stgs.time_limit_secs = 1e-10;
  if (scs_update_settings(w, &stgs) != 0) {
    fail = "update_settings_time_limit: setting the time limit refused";
    goto out;
  }
  exitflag = scs_solve(w, sol, &info, 0);
  if (!upd_hit_time_limit(exitflag, &info)) {
    fail = "update_settings_time_limit: new time limit had no effect";
    goto out;
  }
  if (info.setup_time != setup_time) {
    fail = "update_settings_time_limit: setup time changed";
    goto out;
  }
  /* no call in between: the limit persists */
  exitflag = scs_solve(w, sol, &info, 0);
  if (!upd_hit_time_limit(exitflag, &info)) {
    fail = "update_settings_time_limit: time limit did not persist";
    goto out;
  }

  stgs.time_limit_secs = 0.;
  if (scs_update_settings(w, &stgs) != 0) {
    fail = "update_settings_time_limit: removing the time limit refused";
    goto out;
  }
  exitflag = scs_solve(w, sol, &info, 0);
  if (exitflag != SCS_SOLVED) {
    fail = "update_settings_time_limit: not solved after removing limit";
    goto out;
  }
  scs_finish(w);

  /* a limit passed to scs_init removed before the first solve */
  stgs.time_limit_secs = 1e-10;
  w = scs_init(d, k, &stgs);
  mu_assert("update_settings_time_limit: second scs_init failed", w);
  stgs.time_limit_secs = 0.;
  if (scs_update_settings(w, &stgs) != 0) {
    fail = "update_settings_time_limit: removing the init limit refused";
    goto out;
  }
  exitflag = scs_solve(w, sol, &info, 0);
  if (exitflag != SCS_SOLVED) {
    fail = "update_settings_time_limit: init limit still in force";
    goto out;
  }

out:
  scs_finish(w);
  upd_free(d, k, opt_sol, sol);
  return fail;
}

/* Largest absolute difference between two solutions, over x, y and s. */
static scs_float upd_sol_diff(const ScsSolution *a, const ScsSolution *b,
                              scs_int n, scs_int m) {
  scs_int i;
  scs_float diff = 0.;
  for (i = 0; i < n; ++i) {
    diff = MAX(diff, ABS(a->x[i] - b->x[i]));
  }
  for (i = 0; i < m; ++i) {
    diff = MAX(diff, ABS(a->y[i] - b->y[i]));
    diff = MAX(diff, ABS(a->s[i] - b->s[i]));
  }
  return diff;
}

/*
 * A workspace set up with one set of settings and changed with
 * scs_update_settings before its first solve must solve with the new
 * settings: the stored settings are the new values, the exit flag is that
 * of a workspace set up with them, and (in double precision) the tighter
 * tolerances take more iterations than the base settings. Neither the
 * solutions nor the iteration counts are compared with that workspace's:
 * with USE_LAPACK every vector operation goes through BLAS, and MKL's BLAS
 * can round differently in a third solve than in two identical earlier ones
 * (it did in CI: 275 iterations each, solutions 3.5e-6 apart), and Apple
 * Accelerate is not reproducible run to run (about 1e-5 of the largest entry
 * in double precision; in single precision 425 against 500 iterations on the
 * macOS CI runner). The iterations and the difference are printed. Then
 * max_iters is lowered on the same workspace and the next solve stops
 * there.
 */
static const char *test_update_settings_matches_init(void) {
  ScsData *d;
  ScsCone *k;
  ScsSolution *opt_sol;
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSolution *ref = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSettings base, upd;
  ScsInfo info = {0}, ref_info = {0}, base_info = {0};
  ScsWork *w;
  scs_int exitflag, ref_exitflag;
  scs_float diff, ref_size;
  const char *fail = 0;

  upd_gen_lp(&d, &k, &opt_sol);
  scs_set_default_settings(&base);
  base.verbose = 0;
  upd = base;
  upd.eps_abs = 1e-6;
  upd.eps_rel = 1e-6;
  upd.eps_infeas = 1e-8;
  upd.alpha = 1.8;
  upd.max_iters = 2000;

  ref_exitflag = scs(d, k, &upd, ref, &ref_info);
  scs(d, k, &base, sol, &base_info);

  w = scs_init(d, k, &base);
  mu_assert("update_settings_matches_init: scs_init failed", w);
  if (scs_update_settings(w, &upd) != 0) {
    fail = "update_settings_matches_init: update refused";
    goto out;
  }
  if (w->stgs->eps_abs != upd.eps_abs || w->stgs->eps_rel != upd.eps_rel ||
      w->stgs->eps_infeas != upd.eps_infeas || w->stgs->alpha != upd.alpha ||
      w->stgs->max_iters != upd.max_iters) {
    fail = "update_settings_matches_init: settings not stored";
    goto out;
  }
  exitflag = scs_solve(w, sol, &info, 0);
  diff = upd_sol_diff(sol, ref, d->n, d->m);
  ref_size = MAX(MAX(SCS(norm_inf)(ref->x, d->n), SCS(norm_inf)(ref->y, d->m)),
                 SCS(norm_inf)(ref->s, d->m));
  scs_printf("update_settings_matches_init: %li iterations, %li with the "
             "settings given to scs_init, %li with the base settings; "
             "difference %.2e, largest entry %.2e\n",
             (long)info.iter, (long)ref_info.iter, (long)base_info.iter, diff,
             ref_size);
  if (exitflag != ref_exitflag) {
    fail = "update_settings_matches_init: exit flag differs from scs_init";
    goto out;
  }
#ifndef SFLOAT
  /* In single precision the counts of the two tolerances overlap (with
   * Accelerate the base settings took 350 to 500 iterations and the tight
   * ones 375 to 525), so only double precision checks it, where they are 275
   * and 175 on every backend; the stored settings above and the max_iters
   * update below show the call took effect in both. */
  if (info.iter <= base_info.iter) {
    fail = "update_settings_matches_init: tighter tolerances had no effect";
    goto out;
  }
#endif

  upd.max_iters = 3;
  if (scs_update_settings(w, &upd) != 0) {
    fail = "update_settings_matches_init: max_iters update refused";
    goto out;
  }
  exitflag = scs_solve(w, sol, &info, 0);
  if (info.iter > 3 || exitflag == SCS_SOLVED) {
    fail = "update_settings_matches_init: max_iters update had no effect";
    goto out;
  }

out:
  scs_finish(w);
  SCS(free_sol)(ref);
  upd_free(d, k, opt_sol, sol);
  return fail;
}

/* Each request changes max_iters to 1 together with something that must be
 * refused; if any were applied in part, the solve at the end would stop
 * after one iteration. */
static const char *test_update_settings_refused(void) {
  ScsData *d;
  ScsCone *k;
  ScsSolution *opt_sol;
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSettings base, bad;
  ScsInfo info = {0};
  ScsWork *w;
  scs_int i, exitflag;
  scs_float adapted;
  const char *fail = 0;
  /* one case per setup setting */
  const char *names[] = {"normalize",
                         "scale",
                         "adaptive_scale",
                         "adaptive_diag_scale",
                         "rho_x",
                         "acceleration_lookback",
                         "acceleration_interval",
                         "acceleration_type_1",
                         "acceleration_regularization",
                         "acceleration_relaxation",
                         "write_data_filename"};
  const scs_int n_cases = (scs_int)(sizeof(names) / sizeof(names[0]));

  upd_gen_lp(&d, &k, &opt_sol);
  scs_set_default_settings(&base);
  base.verbose = 0;

  w = scs_init(d, k, &base);
  mu_assert("update_settings_refused: scs_init failed", w);

  for (i = 0; i < n_cases; ++i) {
    bad = base;
    bad.max_iters = 1;
    switch (i) {
    case 0: bad.normalize = !base.normalize; break;
    case 1: bad.scale = 2 * base.scale; break;
    case 2: bad.adaptive_scale = !base.adaptive_scale; break;
    case 3: bad.adaptive_diag_scale = !base.adaptive_diag_scale; break;
    case 4: bad.rho_x = 2 * base.rho_x; break;
    case 5: bad.acceleration_lookback = base.acceleration_lookback + 1; break;
    case 6: bad.acceleration_interval = base.acceleration_interval + 1; break;
    case 7: bad.acceleration_type_1 = !base.acceleration_type_1; break;
    case 8:
      bad.acceleration_regularization = 2 * base.acceleration_regularization;
      break;
    case 9: bad.acceleration_relaxation = 0.5; break;
    case 10: bad.write_data_filename = "update_settings_never_written"; break;
    }
    if (scs_update_settings(w, &bad) >= 0) {
      scs_printf("update_settings_refused: change of %s accepted\n", names[i]);
      fail = "update_settings_refused: change of a setup setting accepted";
      goto out;
    }
  }

#if NO_VALIDATE == 0
  /* values scs_init would refuse */
  for (i = 0; i < 4; ++i) {
    bad = base;
    bad.max_iters = 1;
    switch (i) {
    case 0: bad.time_limit_secs = -1.; break;
    case 1: bad.eps_abs = NAN; break;
    case 2: bad.max_iters = 0; break;
    case 3: bad.alpha = 2.; break;
    }
    if (scs_update_settings(w, &bad) >= 0) {
      fail = "update_settings_refused: invalid value accepted";
      goto out;
    }
  }
#endif

  bad = base;
  bad.max_iters = 1;
  if (scs_update_settings(SCS_NULL, &bad) >= 0 ||
      scs_update_settings(w, SCS_NULL) >= 0) {
    fail = "update_settings_refused: missing input accepted";
    goto out;
  }

  exitflag = scs_solve(w, sol, &info, 0);
  if (info.iter <= 1 || exitflag == SCS_FAILED) {
    fail = "update_settings_refused: a refused update changed the workspace";
    goto out;
  }

  /* The settings given to scs_init are accepted back even when the solver
   * has adapted scale since; the adapted value is kept, not reset. The
   * adaptation is set by hand so that the test does not depend on whether
   * the solve above happened to update the scale. */
  adapted = 10. * base.scale;
  w->stgs->scale = adapted;
  if (scs_update_settings(w, &base) != 0) {
    fail = "update_settings_refused: init settings refused after adaptation";
    goto out;
  }
  if (w->stgs->scale != adapted) {
    fail = "update_settings_refused: adapted scale was overwritten";
    goto out;
  }
  scs_finish(w);

  /* adaptive_diag_scale on with adaptive_scale off: setup turns it off, and
   * the settings given to scs_init are still accepted back */
  base.adaptive_scale = 0;
  base.adaptive_diag_scale = 1;
  w = scs_init(d, k, &base);
  mu_assert("update_settings_refused: second scs_init failed", w);
  if (scs_update_settings(w, &base) != 0) {
    fail = "update_settings_refused: init settings refused after setup "
           "turned adaptive_diag_scale off";
    goto out;
  }

out:
  scs_finish(w);
  upd_free(d, k, opt_sol, sol);
  return fail;
}

/*
 * log_csv_filename is read when each solve starts, so it can be switched on
 * and off between solves. SCS keeps its own copy of the name: the caller's
 * buffer is overwritten after the call.
 */
static const char *test_update_settings_csv(void) {
#if NO_READ_WRITE > 0
  return mu_skipped;
#else
  ScsData *d;
  ScsCone *k;
  ScsSolution *opt_sol;
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsSettings stgs;
  ScsInfo info = {0};
  ScsWork *w;
  FILE *f;
  long size = 0;
  char name[] = "test_update_settings.csv";
  const char *fail = 0;

  remove("test_update_settings.csv");
  upd_gen_lp(&d, &k, &opt_sol);
  scs_set_default_settings(&stgs);
  stgs.verbose = 0;
  w = scs_init(d, k, &stgs);
  mu_assert("update_settings_csv: scs_init failed", w);

  stgs.log_csv_filename = name;
  if (scs_update_settings(w, &stgs) != 0) {
    fail = "update_settings_csv: setting log_csv_filename refused";
    goto out;
  }
  memset(name, 'x', strlen(name)); /* SCS must not use the caller's buffer */
  scs_solve(w, sol, &info, 0);
  f = fopen("test_update_settings.csv", "r");
  if (f) {
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fclose(f);
  }
  remove("test_update_settings.csv");
  if (!f || size <= 0) {
    fail = "update_settings_csv: no log written";
    goto out;
  }

  stgs.log_csv_filename = SCS_NULL;
  if (scs_update_settings(w, &stgs) != 0) {
    fail = "update_settings_csv: clearing log_csv_filename refused";
    goto out;
  }
  scs_solve(w, sol, &info, 0);
  f = fopen("test_update_settings.csv", "r");
  if (f) {
    fclose(f);
    remove("test_update_settings.csv");
    fail = "update_settings_csv: log written after it was switched off";
    goto out;
  }

out:
  scs_finish(w);
  upd_free(d, k, opt_sol, sol);
  return fail;
#endif
}
