#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "scs.h"

#define EXP_CONE_INFINITY_VALUE (1E15)

/*
 * Exponential cone projection routines, from:
 *
 * Projection onto the exponential cone: a univariate root-finding problem,
 *    by Henrik A. Friberg, 2021.
 *
 */

static inline scs_int _isfinite(scs_float x) {
  return ABS(x) < EXP_CONE_INFINITY_VALUE;
}

static inline scs_float _clip(scs_float x, scs_float l, scs_float u) {
  return MAX(l, MIN(u, x));
}

/* As defined in Friberg, 2021 (multiplied by positive polynomial) */
static void hfun(const scs_float *v0, scs_float rho, scs_float *f,
                 scs_float *df) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float exprho = exp(rho);
  scs_float expnegrho = exp(-rho);
  /* function value at v0 */
  *f = ((rho - 1) * r0 + s0) * exprho - (r0 - rho * s0) * expnegrho -
       (rho * (rho - 1) + 1) * t0;
  /* gradient of function at v0 */
  *df = (rho * r0 + s0) * exprho + (r0 - (rho - 1) * s0) * expnegrho -
        (2 * rho - 1) * t0;
}

/* Binary search for the root of the hfun function */
static scs_float root_search_binary(const scs_float *v0, scs_float xl,
                                    scs_float xu, scs_float x) {
#if VERBOSITY > 0
  scs_printf("Exp cone: Newton method failed, resorting to binary search.\n");
#endif
  const scs_float EPS = 1e-12; /* expensive so loosen tol */
  const scs_int MAXITER = 40;
  scs_int i;
  scs_float x_plus = x, f, df;
  for (i = 0; i < MAXITER; i++) {
    hfun(v0, x, &f, &df);
    if (f < 0.0) {
      xl = x;
    } else {
      xu = x;
    }
    /* binary search step */
    x_plus = 0.5 * (xl + xu);
    if (ABS(x_plus - x) <= EPS * MAX(1., ABS(x_plus)) || (x_plus == xl) ||
        (x_plus == xu)) {
      break;
    }
    x = x_plus;
  }
  return x_plus;
}

/* Use damped Newton's to find the root of the hfun function */
static scs_float root_search_newton(const scs_float *v0, scs_float xl,
                                    scs_float xu, scs_float x) {
  /* params taken from Friberg code */
  const scs_float EPS = 1e-15;
  const scs_float DFTOL = 1e-13; /* pow(EPS, 6.0 / 7.0) */
  const scs_int MAXITER = 20;
  const scs_float LODAMP = 0.05;
  const scs_float HIDAMP = 0.95;

  scs_float x_plus, f, df;
  scs_int i;

  for (i = 0; i < MAXITER; i++) {
    hfun(v0, x, &f, &df);

    if (ABS(f) <= EPS) { /* root found */
      break;
    }

    if (f < 0.0) {
      xl = x;
    } else {
      xu = x;
    }

    if (xu <= xl) {
      xu = 0.5 * (xu + xl);
      xl = xu;
      break;
    }

    if (!_isfinite(f) || df < DFTOL) {
      break;
    }

    /* Newton step */
    x_plus = x - f / df;

    if (ABS(x_plus - x) <= EPS * MAX(1., ABS(x_plus))) {
      break;
    }

    if (x_plus >= xu) {
      x = MIN(LODAMP * x + HIDAMP * xu, xu);
    } else if (x_plus <= xl) {
      x = MAX(LODAMP * x + HIDAMP * xl, xl);
    } else {
      x = x_plus;
    }
  }
  if (i < MAXITER) { /* Newton's method converged */
#if VERBOSITY > 0
    scs_printf("Exp cone: Newton iters:%i, f:%.4e, df:%.4e\n", (int)i, f, df);
#endif
    return _clip(x, xl, xu);
  }
  /* Fall back to binary search if Newton failed */
  return root_search_binary(v0, xl, xu, x);
}

/* try heuristic (cheap) projection */
static scs_float proj_primal_exp_cone_heuristic(const scs_float *v0,
                                                scs_float *vp) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float dist, tp, newdist;
  /* perspective boundary */
  vp[2] = MAX(t0, 0);
  vp[1] = 0.0;
  vp[0] = MIN(r0, 0);
  dist = SCS(norm_diff)(v0, vp, 3);

  /* perspective interior */
  if (s0 > 0.0) {
    tp = MAX(t0, s0 * exp(r0 / s0));
    newdist = tp - t0;
    if (newdist < dist) {
      vp[2] = tp;
      vp[1] = s0;
      vp[0] = r0;
      dist = newdist;
    }
  }
  return dist;
}

/* try heuristic (cheap) projection */
static scs_float proj_polar_exp_cone_heuristic(const scs_float *v0,
                                               scs_float *vd) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float dist, td, newdist;
  /* perspective boundary */
  vd[2] = MIN(t0, 0);
  vd[1] = MIN(s0, 0);
  vd[0] = 0.0;
  dist = SCS(norm_diff)(v0, vd, 3);

  /* perspective interior */
  if (r0 > 0.0) {
    td = MIN(t0, -r0 * exp(s0 / r0 - 1));
    newdist = t0 - td;
    if (newdist < dist) {
      vd[2] = td;
      vd[1] = s0;
      vd[0] = r0;
      dist = newdist;
    }
  }
  return dist;
}

static scs_float ppsi(const scs_float *v0) {
  scs_float s0 = v0[1], r0 = v0[0];
  scs_float psi;

  if (r0 > s0) {
    psi = (r0 - s0 + sqrt(r0 * r0 + s0 * s0 - r0 * s0)) / r0;
  } else {
    psi = -s0 / (r0 - s0 - sqrt(r0 * r0 + s0 * s0 - r0 * s0));
  }

  return ((psi - 1) * r0 + s0) / (psi * (psi - 1) + 1);
}

static scs_float pomega(scs_float rho) {
  scs_float val = exp(rho) / (rho * (rho - 1) + 1);

  if (rho < 2.0) {
    val = MIN(val, exp(2.0) / 3);
  }

  return val;
}

static scs_float dpsi(const scs_float *v0) {
  scs_float s0 = v0[1], r0 = v0[0];
  scs_float psi;

  if (s0 > r0) {
    psi = (r0 - sqrt(r0 * r0 + s0 * s0 - r0 * s0)) / s0;
  } else {
    psi = (r0 - s0) / (r0 + sqrt(r0 * r0 + s0 * s0 - r0 * s0));
  }

  return (r0 - psi * s0) / (psi * (psi - 1) + 1);
}

static scs_float domega(scs_float rho) {
  scs_float val = -exp(-rho) / (rho * (rho - 1) + 1);

  if (rho > -1.0) {
    val = MAX(val, -exp(1.0) / 3);
  }

  return val;
}

/* Generate upper and lower search bounds for root of hfun */
static void exp_search_bracket(const scs_float *v0, scs_float pdist,
                               scs_float ddist, scs_float *low_out,
                               scs_float *upr_out) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float baselow = -EXP_CONE_INFINITY_VALUE,
            baseupr = EXP_CONE_INFINITY_VALUE;
  scs_float low = -EXP_CONE_INFINITY_VALUE, upr = EXP_CONE_INFINITY_VALUE;

  scs_float Dp = SQRTF(pdist * pdist - MIN(s0, 0) * MIN(s0, 0));
  scs_float Dd = SQRTF(ddist * ddist - MIN(r0, 0) * MIN(r0, 0));

  scs_float curbnd, fl, fu, df, tpu, tdl;

  if (t0 > 0) {
    curbnd = log(t0 / ppsi(v0));
    low = MAX(low, curbnd);
  } else if (t0 < 0) {
    curbnd = -log(-t0 / dpsi(v0));
    upr = MIN(upr, curbnd);
  }

  if (r0 > 0) {
    baselow = 1 - s0 / r0;
    low = MAX(low, baselow);
    tpu = MAX(1e-12, MIN(Dd, Dp + t0));
    curbnd = MAX(low, baselow + tpu / r0 / pomega(low));
    upr = MIN(upr, curbnd);
  }

  if (s0 > 0) {
    baseupr = r0 / s0;
    upr = MIN(upr, baseupr);
    tdl = -MAX(1e-12, MIN(Dp, Dd - t0));
    curbnd = MIN(upr, baseupr - tdl / s0 / domega(upr));
    low = MAX(low, curbnd);
  }

  /* Guarantee valid bracket */
  /* TODO do we need these 2 lines? */
  low = _clip(MIN(low, upr), baselow, baseupr);
  upr = _clip(MAX(low, upr), baselow, baseupr);

  if (low != upr) {
    hfun(v0, low, &fl, &df);
    hfun(v0, upr, &fu, &df);

    if (fl * fu > 0) {
      if (ABS(fl) < ABS(fu)) {
        upr = low;
      } else {
        low = upr;
      }
    }
  }

  *low_out = low;
  *upr_out = upr;
}

/* convert from rho to primal projection */
static scs_float proj_sol_primal_exp_cone(const scs_float *v0, scs_float rho,
                                          scs_float *vp) {
  scs_float linrho = (rho - 1) * v0[0] + v0[1];
  scs_float exprho = exp(rho);
  scs_float quadrho, dist;
  if (linrho > 0 && _isfinite(exprho)) {
    quadrho = rho * (rho - 1) + 1;
    vp[2] = exprho * linrho / quadrho;
    vp[1] = linrho / quadrho;
    vp[0] = rho * linrho / quadrho;
    dist = SCS(norm_diff)(vp, v0, 3);
  } else {
    vp[2] = EXP_CONE_INFINITY_VALUE;
    vp[1] = 0.0;
    vp[0] = 0.0;
    dist = EXP_CONE_INFINITY_VALUE;
  }
  return dist;
}

/* convert from rho to polar projection */
static scs_float proj_sol_polar_exp_cone(const scs_float *v0, scs_float rho,
                                         scs_float *vd) {
  scs_float linrho = v0[0] - rho * v0[1];
  scs_float exprho = exp(-rho);
  scs_float quadrho, dist;
  if (linrho > 0 && _isfinite(exprho)) {
    quadrho = rho * (rho - 1) + 1;
    vd[2] = -exprho * linrho / quadrho;
    vd[1] = (1 - rho) * linrho / quadrho;
    vd[0] = linrho / quadrho;
    dist = SCS(norm_diff)(v0, vd, 3);
  } else {
    vd[2] = -EXP_CONE_INFINITY_VALUE;
    vd[1] = 0.0;
    vd[0] = 0.0;
    dist = EXP_CONE_INFINITY_VALUE;
  }
  return dist;
}

static inline void _copy(scs_float *dest, scs_float *src) {
  dest[0] = src[0];
  dest[1] = src[1];
  dest[2] = src[2];
}

/* Project onto primal or dual exponential cone, performed in-place.
 * If `primal=0` then project on the dual cone, otherwise project
 * onto primal. Taken from algorithm in Friberg, 2021.
 */
scs_float SCS(proj_pd_exp_cone)(scs_float *v0, scs_int primal) {
  scs_float TOL = 1e-8; /* pow(1e-10, 2.0 / 3.0); */
  scs_float xl, xh, pdist, ddist, err, rho, dist_hat;
  scs_float vp[3], vd[3], v_hat[3];
  scs_int opt;
  if (!primal) {
    /* This routine actually projects onto primal and polar cones
     * simultaneously. So to make it project onto dual, use this:
     * Pi_{C^*}(v0) = -Pi_{C^polar}(-v0)
     */
    v0[0] *= -1.;
    v0[1] *= -1.;
    v0[2] *= -1.;
  }

  pdist = proj_primal_exp_cone_heuristic(v0, vp);
  ddist = proj_polar_exp_cone_heuristic(v0, vd);

  err = ABS(vp[0] + vd[0] - v0[0]);
  err = MAX(err, ABS(vp[1] + vd[1] - v0[1]));
  err = MAX(err, ABS(vp[2] + vd[2] - v0[2]));

  /* Skip root search if presolve rules apply
   * or optimality conditions are satisfied
   */
  opt = (v0[1] <= 0 && v0[0] <= 0);
  opt |= (MIN(pdist, ddist) <= TOL);
  opt |= (err <= TOL && SCS(dot)(vp, vd, 3) <= TOL);
  if (opt) {
    if (primal) {
      _copy(v0, vp);
      return pdist;
    }
    /* polar cone -> dual cone */
    v0[0] = -vd[0];
    v0[1] = -vd[1];
    v0[2] = -vd[2];
    return ddist;
  }

  exp_search_bracket(v0, pdist, ddist, &xl, &xh);
  rho = root_search_newton(v0, xl, xh, 0.5 * (xl + xh));

  if (primal) {
    /* primal cone projection */
    dist_hat = proj_sol_primal_exp_cone(v0, rho, v_hat);
    if (dist_hat <= pdist) {
      _copy(vp, v_hat);
      pdist = dist_hat;
    }
    _copy(v0, vp);
    return pdist;
  }
  /* polar cone projection */
  dist_hat = proj_sol_polar_exp_cone(v0, rho, v_hat);
  if (dist_hat <= ddist) {
    _copy(vd, v_hat);
    ddist = dist_hat;
  }
  /* polar cone -> dual cone */
  v0[0] = -vd[0];
  v0[1] = -vd[1];
  v0[2] = -vd[2];
  return ddist;
}
