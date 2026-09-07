#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

/*
 * Spectral cones under mid-solve metric changes.
 *
 * Unlike the closed-form projections of the standard cones, the spectral-cone
 * projections are iterative inner solvers that carry warm-start state across
 * outer iterations (log_cone_warmstarts, saved_log_projs). Both the scalar
 * adaptive_scale and the dynamic diagonal rescaling change the diag_r metric
 * mid-solve, and that state belongs to the metric that produced it, so
 * update_scale drops it via SCS(reset_cone_cache). This test covers the
 * interaction that hook exists for.
 *
 * A regression here is silent: the projections are the only writers of the
 * dual iterate and the SCS residuals never test cone membership, so a bad
 * projection terminates 'solved' with a correct objective and a y far outside
 * the dual cone. The objective checks below would not catch it; the cone
 * distance checks are the point of this test.
 *
 * The problem is the exp_design D-optimal design instance (one logdet cone,
 * the projection with the Newton/IPM inner solver and warm-start cache). We
 * solve it repeatedly from initial scales spanning eight orders of magnitude,
 * which forces many metric updates rather than the handful a well-scaled
 * start produces, and assert that every solve lands on the cone.
 */
static const char *test_spectral_metric_updates(void) {
  /* initial scales far from the value the heuristic settles on, so that
   * adaptive_scale and the diagonal rescaling both fire repeatedly.
   * The range stops at 1e-3: this instance's logdet projection diverges
   * from a 1e-4 start (on master too, at iteration 0 before any metric
   * update), which is a separate robustness limit and not what this test
   * is about. */
  static const scs_float scales[] = {1e-3, 1e-2, 1e-1, 1e0,
                                     1e1,  1e2,  1e3,  1e4};
  const scs_int n_scales = (scs_int)(sizeof(scales) / sizeof(scales[0]));
  scs_int total_scale_updates = 0;
  scs_int trial;

  /* exp_design data: m = 15, n = 8, one logdet cone of dimension 3 */
  scs_float Ax[] = {
      -1.,         1.,          -1.,         3.24,        0.16,
      1.,          4.84,        3.61,        1.,          -1.,
      2.54558441,  -0.11313708, -0.14142136, 1.24450793,  0.26870058,
      -2.12132034, -1.,         2.03646753,  0.05656854,  0.56568542,
      0.93338095,  4.03050865,  0.28284271,  -1.,         1.,
      0.04,        0.01,        0.16,        0.01,        2.25,
      -1.,         1.13137085,  -0.02828427, -0.05656854, 0.16970563,
      0.21213203,  -0.42426407, -1.,         0.64,        0.01,
      0.16,        0.09,        2.25,        0.04,        -1.};
  scs_int Ai[] = {7,  0,  8, 1, 2, 3, 4, 5,  6,  9, 1, 2, 3, 4, 5,
                  6,  10, 1, 2, 3, 4, 5, 6,  11, 1, 2, 3, 4, 5, 6,
                  12, 1,  2, 3, 4, 5, 6, 13, 1,  2, 3, 4, 5, 6, 14};
  scs_int Ap[] = {0, 1, 3, 10, 17, 24, 31, 38, 45};

  scs_float b[] = {1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.};
  scs_float c[] = {1., 0., 0., 0., 0., 0., 0., 0.};

  scs_int m = 15;
  scs_int n = 8;
  scs_int d_array[] = {3};

  /* computed using mosek (Ax above is truncated, mosek solved the problem
   * with the non-truncated data) */
  scs_float opt = 3.0333290743428574;

  for (trial = 0; trial < n_scales; ++trial) {
    ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
    ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
    ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
    ScsConeWork *cone_work;
    ScsInfo info = {0};
    scs_int exitflag;
    scs_float perr, derr, ydist, sdist;

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

    k->z = 1;
    k->l = 6;
    k->d = d_array;
    k->dsize = 1;

    scs_set_default_settings(stgs);
    stgs->eps_abs = 1e-7;
    stgs->eps_rel = 1e-7;
    stgs->eps_infeas = 1e-9;
    stgs->scale = scales[trial];
    /* both metric-update mechanisms on: the scalar heuristic and the
     * dynamic diagonal rescaling that updates far more often */
    stgs->adaptive_scale = 1;
    stgs->adaptive_diag_scale = 1;

    exitflag = scs(d, k, stgs, sol, &info);

    perr = SCS(dot)(d->c, sol->x, d->n) - opt;
    derr = -SCS(dot)(d->b, sol->y, d->m) - opt;

    scs_printf("scale %.0e: %li scale updates, primal err %.4e, dual err "
               "%.4e\n",
               stgs->scale, (long)info.scale_updates, perr, derr);

    mu_assert("test_spectral_metric_updates: SCS failed to produce "
              "outputflag SCS_SOLVED",
              exitflag == SCS_SOLVED);
    mu_assert_less("test_spectral_metric_updates: primal obj ERROR", ABS(perr),
                   1e-4);
    mu_assert_less("test_spectral_metric_updates: dual obj ERROR", ABS(derr),
                   1e-4);

    /* the actual regression check: a projection that mishandled the metric
     * change leaves y outside the dual cone (and s outside the cone) while
     * the objective above still looks right */
    cone_work = SCS(init_cone)(k, m);
    mu_assert("test_spectral_metric_updates: failed to init cone work",
              cone_work != SCS_NULL);
    ydist = get_dual_cone_dist(sol->y, cone_work, m);
    sdist = get_pri_cone_dist(sol->s, cone_work, m);
    SCS(finish_cone)(cone_work);

    scs_printf("scale %.0e: y cone dist %.4e, s cone dist %.4e\n",
               stgs->scale, ydist, sdist);

    mu_assert_less("test_spectral_metric_updates: y cone dist ERROR",
                   ABS(ydist), 1e-5);
    mu_assert_less("test_spectral_metric_updates: s cone dist ERROR",
                   ABS(sdist), 1e-5);

    total_scale_updates += info.scale_updates;

    scs_free(d->A);
    scs_free(k);
    scs_free(stgs);
    scs_free(d);
    SCS(free_sol)(sol);
  }

  /* Guard the guard: if the metric never changed mid-solve, the checks above
   * are vacuous and this test has silently stopped covering the interaction.
   */
  mu_assert("test_spectral_metric_updates: no metric updates were triggered, "
            "so the test no longer exercises mid-solve metric changes",
            total_scale_updates > 0);

  return 0;
}
