#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * The box cone projection (proj_box_cone in cones.c), called through
 * SCS(proj_dual_cone) as SCS calls it. It projects (t0, x) onto
 *   K = { (t, s) : t l <= s <= t u, t >= 0 }
 * in the metric rho_t (t - t0)^2 + sum_j r_j (s_j - x_j)^2: t* minimises a
 * convex piecewise quadratic whose derivative is
 *   g(t) = rho_t (t - t0) + sum_{x_j > t u_j} r_j (t u_j - x_j) u_j
 *                         + sum_{x_j < t l_j} r_j (t l_j - x_j) l_j,
 * continuous and increasing, and s = x clipped to [t l, t u].
 *
 * proj_dual_cone(x) returns x + P(-x), P the projection onto K in the metric
 * diag(1 / r_y) (rho_t = 1 / r_y[0], r_j = 1 / r_y[j + 1]) of the input
 * -r_y x, so the projection of a point p is recovered from the call on
 * x = -p / r_y as r_y (returned - x). The cone work keeps the last t as the
 * next call's warm start, as in a solve.
 */

/* g(t) for the projection of (t0, x) (see above); rinv[j] = r_j */
static scs_float box_proj_g(scs_float t, scs_float t0, const scs_float *x,
                            const scs_float *bl, const scs_float *bu,
                            scs_int nb, scs_float rho_t,
                            const scs_float *rinv) {
  scs_float g = rho_t * (t - t0);
  scs_int j;
  for (j = 0; j < nb; ++j) {
    if (x[j] > t * bu[j]) {
      g += rinv[j] * (t * bu[j] - x[j]) * bu[j];
    } else if (x[j] < t * bl[j]) {
      g += rinv[j] * (t * bl[j] - x[j]) * bl[j];
    }
  }
  return g;
}

/* Project p = (t0, x) (length nb + 1) with weights r_y in place through
 * proj_dual_cone, as described above. The round trip through x + P(-x)
 * costs about eps max(|t0|, |t|) of absolute accuracy in t. */
static scs_int box_proj_call(ScsConeWork *cw, scs_float *p,
                             const scs_float *r_y, scs_float *work,
                             scs_int nb) {
  scs_int i, ret;
  for (i = 0; i <= nb; ++i) {
    work[i] = -p[i] / r_y[i];
  }
  ret = SCS(proj_dual_cone)(work, cw, SCS_NULL, (scs_float *)r_y);
  for (i = 0; i <= nb; ++i) {
    p[i] = r_y[i] * (work[i] + p[i] / r_y[i]);
  }
  return ret;
}

/* SCS 3.3's Newton loop stopped when a step fell below 1e-11 max(t, 1). From
 * a warm start at 0, a row with a large bound (here 1e12; qafiro_tiny_qp
 * without normalization has bounds of 1e20) and x_j > 0 is clipped at t = 0,
 * and its curvature r u^2 = 1e24 makes the first step x_j / u_j = 1e-12,
 * which met that test: the loop returned t = 1e-12 where t* = 0.5, and with
 * the next call warm-started there it could not leave. A projection that
 * stops short like that returns t off by about t*; the test allows 1e-6. */
static const char *test_box_proj_kink_trap(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsConeWork *cw;
  scs_float bl[] = {0., -1.}, bu[] = {1e12, 1.};
  scs_float r_y[] = {1., 1., 1.}, rinv[] = {1., 1.};
  /* t* = 0 for the first input, which leaves the warm start at 0 */
  scs_float p0[] = {-1., 0., 0.};
  /* t* = 0.5: rows 1 and 2 are free for t > 1e-12, so g(t) = t - 0.5 */
  scs_float p1[] = {0.5, 1., 0.};
  scs_float work[3], t_star = 0.5;
  const char *fail = SCS_NULL;

  k->bsize = 3;
  k->bl = bl;
  k->bu = bu;
  cw = SCS(init_cone)(k, 3);
  mu_assert("test_box_proj_kink_trap: cone init failed", cw);
  box_proj_call(cw, p0, r_y, work, 2);
  if (ABS(p0[0]) > 1e-12) {
    fail = "test_box_proj_kink_trap: t of the first projection is not 0";
  }
  box_proj_call(cw, p1, r_y, work, 2);
  scs_printf("test_box_proj_kink_trap: t = %.6e (t* = %.1f), g(t) = %.3e\n",
             p1[0], t_star,
             box_proj_g(p1[0], 0.5, &(p1[1]), bl, bu, 2, 1., rinv));
  SCS(finish_cone)(cw);
  scs_free(k);
  if (fail) {
    return fail;
  }
  mu_assert_less("test_box_proj_kink_trap: |t - t*|", ABS(p1[0] - t_star),
                 1e-6);
  return 0;
}

/* deterministic pseudo-random numbers in [0, 1) (64-bit LCG, Knuth's
 * MMIX constants), so the instances are the same on every platform */
static scs_float box_proj_rand(unsigned long long *state) {
  *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
  return (scs_float)((double)(*state >> 11) / 9007199254740992.);
}

#define BOX_PROJ_CONES (300)
#define BOX_PROJ_CALLS (5)
#define BOX_PROJ_MAX_ROWS (40)

/* Optimality of the returned point on 300 cone works of up to 40 rows, five
 * calls each (warm-started from the previous call, as in a solve), with
 * random weights and bounds: finite, infinite, zero width, and steep (bounds
 * of 1e6 to 1e14 with their kinks near 0, the 3.3 loop's trap; larger ones
 * would overflow u^2 in single precision), and inputs with kinks within a
 * few ulps of t*. The returned t is t* to within d
 * max(t, 1) if g changes sign across [t - d max(t, 1), t + d max(t, 1)]
 * (only g(t + d max(t, 1)) > 0 when that interval reaches below 0), and s
 * must be x clipped to [t' l, t' u] for a t' in that interval, to d max(|x|,
 * 1). d = 1e-7 (1e-3 in single precision) is far
 * above the rounding of these well-scaled instances and far below what a
 * projection stopped short of its root gives (the 3.3 loop: t off by up to
 * 100% of t* on the steep instances). */
static const char *test_box_proj_optimality(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  scs_float bl[BOX_PROJ_MAX_ROWS], bu[BOX_PROJ_MAX_ROWS];
  scs_float r_y[BOX_PROJ_MAX_ROWS + 1], rinv[BOX_PROJ_MAX_ROWS];
  scs_float p[BOX_PROJ_MAX_ROWS + 1], x[BOX_PROJ_MAX_ROWS];
  scs_float work[BOX_PROJ_MAX_ROWS + 1];
  scs_float d = VERIFY_TOL(1e-7, 1e-3), t, t0, h, lo, hi, s, e;
  scs_float worst = 0., one = 1.;
  unsigned long long state = 20261119ULL;
  scs_int i, j, call, nb, family, nbad = 0, ncalls = 0;
  const char *fail = SCS_NULL;

  for (i = 0; i < BOX_PROJ_CONES && !fail; ++i) {
    ScsConeWork *cw;
    nb = 1 + (scs_int)(box_proj_rand(&state) * BOX_PROJ_MAX_ROWS);
    family = i % 3; /* 0 generic, 1 steep rows, 2 kinks at the root */
    for (j = 0; j < nb; ++j) {
      scs_float a = box_proj_rand(&state), b = box_proj_rand(&state);
      bl[j] = a < 0.15 ? -INFINITY : (a < 0.3 ? 0. : -5. * b);
      a = box_proj_rand(&state);
      b = box_proj_rand(&state);
      bu[j] = a < 0.15 ? INFINITY : (a < 0.25 ? bl[j] : 5. * b);
      if (bu[j] < bl[j]) { /* keep l <= u */
        bu[j] = bl[j] + 5. * b;
      }
      if (family == 1 && box_proj_rand(&state) < 0.3) {
        /* steep: a bound of 1e6 to 1e14, its kink x / u near 0 */
        bu[j] = POWF(10., 6. + 8. * box_proj_rand(&state));
        bl[j] = box_proj_rand(&state) < 0.5 ? -INFINITY : 0.;
      }
    }
    k->bsize = nb + 1;
    k->bl = bl;
    k->bu = bu;
    cw = SCS(init_cone)(k, nb + 1);
    if (!cw) {
      fail = "test_box_proj_optimality: cone init failed";
      break;
    }
    for (call = 0; call < BOX_PROJ_CALLS; ++call) {
      /* weights: rho_t = 1 / r_y[0] from 1e-3 to 1e3, r_j from 1e-2 to 1e2 */
      r_y[0] = POWF(10., -3. + 6. * box_proj_rand(&state));
      for (j = 0; j < nb; ++j) {
        r_y[j + 1] = POWF(10., -2. + 4. * box_proj_rand(&state));
        rinv[j] = 1. / r_y[j + 1];
      }
      t0 = -2. + 5. * box_proj_rand(&state);
      for (j = 0; j < nb; ++j) {
        x[j] = -10. + 20. * box_proj_rand(&state);
        if (family == 1 && bu[j] > 1e5) {
          x[j] = 10. * box_proj_rand(&state);
        }
      }
      if (family == 2) {
        /* t* = t0 up to a few ulps: every row is free at t0, or has its kink
         * x / u within a few ulps of t0 (its term in g is then 0 up to that
         * rounding) */
        t0 = 0.1 + 2. * box_proj_rand(&state);
        for (j = 0; j < nb; ++j) {
          lo = bl[j] > -INFINITY ? t0 * bl[j] : t0 * MIN(bu[j], 0.) - 10.;
          hi = bu[j] < INFINITY ? t0 * bu[j] : t0 * MAX(bl[j], 0.) + 10.;
          x[j] = lo + (hi - lo) * box_proj_rand(&state);
          if (bu[j] > 0. && bu[j] < INFINITY && box_proj_rand(&state) < 0.5) {
            x[j] = t0 * bu[j] * (1. + (box_proj_rand(&state) - 0.5) * 1e-14);
          }
        }
      }
      p[0] = t0;
      for (j = 0; j < nb; ++j) {
        p[j + 1] = x[j];
      }
      box_proj_call(cw, p, r_y, work, nb);
      ncalls++;
      t = p[0];
      h = d * MAX(t, one);
      lo = t - h;
      hi = t + h;
      e = 0.;
      /* g(hi) > 0, and g(lo) < 0 unless lo < 0 */
      if (!(box_proj_g(hi, t0, x, bl, bu, nb, 1. / r_y[0], rinv) > 0.) ||
          (lo > 0. && !(box_proj_g(lo, t0, x, bl, bu, nb, 1. / r_y[0],
                                   rinv) < 0.))) {
        e = 1.;
      }
      for (j = 0; j < nb; ++j) {
        /* x clipped to [t l, t u] (an infinite bound never clips, also at
         * t = 0); s may be that clip at any t within h of the returned one,
         * which moves it by at most h |bound| */
        s = MIN(MAX(x[j], bl[j] > -INFINITY ? t * bl[j] : -INFINITY),
                bu[j] < INFINITY ? t * bu[j] : INFINITY);
        lo = MAX(bl[j] > -INFINITY ? ABS(bl[j]) : 0.,
                 bu[j] < INFINITY ? ABS(bu[j]) : 0.);
        e = MAX(e, (ABS(p[j + 1] - s) - h * lo) / MAX(ABS(x[j]), one));
      }
      worst = MAX(worst, e);
      if (e > d) {
        if (nbad++ < 5) {
          scs_printf("test_box_proj_optimality: cone %li (family %li, %li "
                     "rows) call %li: t = %.17g, error %.3e\n",
                     (long)i, (long)family, (long)nb, (long)call, t, e);
        }
      }
    }
    SCS(finish_cone)(cw);
  }
  scs_free(k);
  scs_printf("test_box_proj_optimality: %li projections, %li not optimal, "
             "largest error %.3e\n",
             (long)ncalls, (long)nbad, worst);
  if (fail) {
    return fail;
  }
  mu_assert("test_box_proj_optimality: a projection is not optimal",
            nbad == 0);
  return 0;
}
