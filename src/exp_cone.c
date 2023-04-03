#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "scs.h"
#include <string.h>

#define EXP_CONE_INF (1E15)

/*
 * Exponential cone projection routines, from:
 *
 * Projection onto the exponential cone: a univariate root-finding problem,
 *    by Henrik A. Friberg, 2021.
 *
 */

static inline scs_int _isfinite(float x) {
  return ABS(x) < EXP_CONE_INF;
}

static inline float _clip(float x, float l, float u) {
  return MIN(MAX(x, l), u);
}

static void hfun(const scs_float *v0, scs_float rho, scs_float *f, scs_float *df) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float exprho = exp(rho);
  scs_float expnegrho = exp(-rho);
  *f = ((rho - 1) * r0 + s0) * exprho - (r0 - rho * s0) * expnegrho -
       (rho * (rho - 1) + 1) * t0;
  *df = (rho * r0 + s0) * exprho + (r0 - (rho - 1) * s0) * expnegrho -
        (2 * rho - 1) * t0;
}

static scs_float root_binary_search(const scs_float *farg, scs_float xl,
                                    scs_float xh, scs_float x0) {
  const scs_float EPS = 1e-15;
  scs_int MAXITER = 20, i;
  scs_float xx, f, df;
  for (i = 0; i < MAXITER; i++) {
    hfun(farg, x0, &f, &df);
    if (f < 0.0) {
      xl = x0;
    } else {
      xh = x0;
    }

    xx = 0.5 * (xl + xh);
    if (ABS(xx - x0) <= EPS * MAX(1., ABS(xx)) || (xx == xl) || (xx == xh)) {
      break;
    }
  }
  return xx;
}

static scs_float root_search_newton(const scs_float *farg, scs_float xl,
                                    scs_float xh, scs_float x0) {
  const scs_float EPS = 1e-15;
  const scs_float DFTOL = pow(EPS, 6.0 / 7.0);
  const scs_int MAXITER = 20;
  const scs_float LODAMP = 0.05;
  const scs_float HIDAMP = 0.95;

  scs_float xx = x0;
  scs_float f, df;

  scs_int i;

  for (i = 0; i < MAXITER; i++) {
    hfun(farg, x0, &f, &df);
    if (f < 0.0) {
      xl = x0;
    } else {
      xh = x0;
    }

    if (xh <= xl) {
      break;
    }

    if (_isfinite(f) && df >= DFTOL) {
      xx = x0 - f / df;
    } else {
      break;
    }

    if (ABS(xx - x0) <= EPS * MAX(1., ABS(xx))) {
      break;
    }

    if (xx >= xh) {
      x0 = MIN(LODAMP * x0 + HIDAMP * xh, xh);
    } else if (xx <= xl) {
      x0 = MAX(LODAMP * x0 + HIDAMP * xl, xl);
    } else {
      x0 = xx;
    }
  }

  if (i < MAXITER) { /* newton method converged */
    return MAX(xl, MIN(xh, xx));
  }
  return root_binary_search(farg, xl, xh, x0);
}

static scs_float proj_primal_exp_cone_heuristic(const scs_float *v0,
                                                scs_float *vp) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float dist;
  /* perspective boundary */
  vp[2] = MAX(t0, 0);
  vp[1] = 0.0;
  vp[0] = MIN(r0, 0);
  dist = SCS(norm_diff)(v0, vp, 3);

  /* perspective scs_interior */
  if (s0 > 0.0) {
    scs_float tp = MAX(t0, s0 * exp(r0 / s0));
    scs_float newdist = tp - t0;
    if (newdist < dist) {
      vp[2] = tp;
      vp[1] = s0;
      vp[0] = r0;
      dist = newdist;
    }
  }
  return dist;
}

static scs_float proj_polar_exp_cone_heuristic(const scs_float *v0,
                                               scs_float *vd) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float dist;
  /* perspective boundary */
  vd[2] = MIN(t0, 0);
  vd[1] = MIN(s0, 0);
  vd[0] = 0.0;
  dist = SCS(norm_diff)(v0, vd, 3);

  /* perspective scs_interior */
  if (r0 > 0.0) {
    scs_float td = MIN(t0, -r0 * exp(s0 / r0 - 1));
    scs_float newdist = t0 - td;
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

static void exp_search_bracket(const scs_float *v0, scs_float pdist,
                               scs_float ddist, scs_float *low_out,
                               scs_float *upr_out) {
  scs_float t0 = v0[2], s0 = v0[1], r0 = v0[0];
  scs_float baselow = -EXP_CONE_INF, baseupr = EXP_CONE_INF;
  scs_float low = -EXP_CONE_INF, upr = EXP_CONE_INF;

  scs_float Dp = sqrt(pdist * pdist - MIN(s0, 0) * MIN(s0, 0));
  scs_float Dd = sqrt(ddist * ddist - MIN(r0, 0) * MIN(r0, 0));

  if (t0 > 0) {
    scs_float tpl = t0;
    scs_float curbnd = log(tpl / ppsi(v0));
    low = MAX(low, curbnd);
  }

  if (t0 < 0) {
    scs_float tdu = t0;
    scs_float curbnd = -log(-tdu / dpsi(v0));
    upr = MIN(upr, curbnd);
  }

  if (r0 > 0) {
    baselow = 1 - s0 / r0;
    low = MAX(low, baselow);

    scs_float tpu = MAX(1e-12, MIN(Dd, Dp + t0));
    scs_float palpha = low;
    scs_float curbnd = MAX(palpha, baselow + tpu / r0 / pomega(palpha));
    upr = MIN(upr, curbnd);
  }

  if (s0 > 0) {
    baseupr = r0 / s0;
    upr = MIN(upr, baseupr);

    scs_float tdl = -MAX(1e-12, MIN(Dp, Dd - t0));
    scs_float dalpha = upr;
    scs_float curbnd = MIN(dalpha, baseupr - tdl / s0 / domega(dalpha));
    low = MAX(low, curbnd);
  }

  /* Guarantee valid bracket */
  /* TODO do we need these 4 lines? */
  low = MIN(low, upr);
  upr = MAX(low, upr);
  low = _clip(low, baselow, baseupr);
  upr = _clip(upr, baselow, baseupr);

  if (low != upr) {
    scs_float fl, fu;
    scs_float df;
    hfun(v0, low, &fl, &df);
    hfun(v0, upr, &fu, &df);

    if (!(fl * fu < 0)) {
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
    vp[2] = EXP_CONE_INF;
    vp[1] = 0.0;
    vp[0] = 0.0;
    dist = EXP_CONE_INF;
  }
  return dist;
}

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
    vd[2] = -EXP_CONE_INF;
    vd[1] = 0.0;
    vd[0] = 0.0;
    dist = EXP_CONE_INF;
  }
  return dist;
}

/* project onto primal or dual exponential conem performed in-place */
scs_float SCS(proj_pd_exp_cone)(scs_float *v0, scs_int primal) {
  scs_float TOL = pow(1e-10, 2.0 / 3.0);
  scs_float xl, xh;
  scs_float vp[3], vd[3];
  if (!primal) {
    v0[0] *= -1.;
    v0[1] *= -1.;
    v0[2] *= -1.;
  }

  scs_float pdist = proj_primal_exp_cone_heuristic(v0, vp);
  scs_float ddist = proj_polar_exp_cone_heuristic(v0, vd);

  scs_float err = ABS(vp[0] + vd[0] - v0[0]);
  err = MAX(err, ABS(vp[1] + vd[1] - v0[1]));
  err = MAX(err, ABS(vp[2] + vd[2] - v0[2]));

  /* Skip root search if presolve rules apply
   * or optimality conditions are satisfied
   */
  scs_int opt = (v0[1] <= 0 && v0[0] <= 0);
  opt |= (MIN(pdist, ddist) <= TOL);
  opt |= (err <= TOL && SCS(dot)(vp, vd, 3) <= TOL);
  if (opt) {
    if (primal) {
      memcpy(v0, vp, 3 * sizeof(scs_float));
      return pdist;
    }
    memcpy(v0, vd, 3 * sizeof(scs_float));
    /* polar cone -> dual cone */
    v0[0] *= -1.;
    v0[1] *= -1.;
    v0[2] *= -1.;
    return ddist;
  }

  exp_search_bracket(v0, pdist, ddist, &xl, &xh);
  scs_float rho = root_search_newton(v0, xl, xh, 0.5 * (xl + xh));

  if (primal) { /* primal cone projection */
    scs_float vp1[3], pdist1;
    pdist1 = proj_sol_primal_exp_cone(v0, rho, vp1);
    if (pdist1 <= pdist) {
      memcpy(vp, vp1, 3 * sizeof(scs_float));
      pdist = pdist1;
    }
    memcpy(v0, vp, 3 * sizeof(scs_float));
    return pdist;
  } /* polar cone projection */
  scs_float vd1[3], ddist1;
  ddist1 = proj_sol_polar_exp_cone(v0, rho, vd1);
  if (ddist1 <= ddist) {
    memcpy(vd, vd1, 3 * sizeof(scs_float));
    ddist = ddist1;
  }
  memcpy(v0, vd, 3 * sizeof(scs_float));
  /* polar cone -> dual cone */
  v0[0] *= -1.;
  v0[1] *= -1.;
  v0[2] *= -1.;
  return ddist;
}

