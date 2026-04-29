#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "scs.h"

/*
 * Test that the stable fused single-loop root_plus computation gives the same
 * result as the original five-separate-dot_r formulation on well-conditioned
 * cases and improves a cancellation-prone quadratic.
 *
 * root_plus solves: tau = (-b + sqrt(b^2 - 4ac)) / (2a)
 *   where a = tau_scale + g'Rg
 *         b = mu'Rg - 2 p'Rg - eta * tau_scale
 *         c = p'Rp - p'Rmu
 *   and x'Ry = sum_i x[i]*y[i]*r[i]  (R-weighted dot product)
 *
 * The old code computed a, b, c via five separate dot_r calls.
 * The new code computes all five weighted dot products in one fused loop.
 * Both must produce the same quadratic formula result.
 */

/* Old implementation: separate weighted dot products */
static scs_float dot_r_ref(const scs_float *x, const scs_float *y,
                           const scs_float *r, scs_int nm) {
  scs_int i;
  scs_float ip = 0.0;
  for (i = 0; i < nm; ++i) {
    ip += x[i] * y[i] * r[i];
  }
  return ip;
}

static scs_float root_plus_old(const scs_float *g, const scs_float *p,
                               const scs_float *mu, const scs_float *r,
                               scs_int nm, scs_float tau_scale,
                               scs_float eta) {
  scs_float a, b, c, rad;
  a = tau_scale + dot_r_ref(g, g, r, nm);
  b = dot_r_ref(mu, g, r, nm) - 2 * dot_r_ref(p, g, r, nm) -
      eta * tau_scale;
  c = dot_r_ref(p, p, r, nm) - dot_r_ref(p, mu, r, nm);
  rad = b * b - 4 * a * c;
  return (-b + SQRTF(MAX(rad, 0.))) / (2 * a);
}

static scs_float root_plus_stable_coeffs(scs_float a, scs_float b,
                                         scs_float c) {
  scs_float rad, sqrt_rad, q;
  if (!isfinite(a) || !isfinite(b) || !isfinite(c) || a <= 0.) {
    return NAN;
  }
  rad = b * b - 4 * a * c;
  if (!isfinite(rad)) {
    return NAN;
  }
  if (rad < 0.) {
    return -b / (2 * a);
  }
  sqrt_rad = SQRTF(rad);
  if (b <= 0.) {
    return (-b + sqrt_rad) / (2 * a);
  }
  q = -0.5 * (b + sqrt_rad);
  return q != 0. ? c / q : 0.;
}

/* New implementation: fused single loop with stable quadratic formula */
static scs_float root_plus_new(const scs_float *g, const scs_float *p,
                               const scs_float *mu, const scs_float *r,
                               scs_int nm, scs_float tau_scale,
                               scs_float eta) {
  scs_int i;
  scs_float gg = 0., mug = 0., pg = 0., pp = 0., pmu = 0.;
  scs_float a, b, c;
  for (i = 0; i < nm; ++i) {
    scs_float ri = r[i], gi = g[i], pi = p[i], mui = mu[i];
    gg  += gi  * gi  * ri;
    mug += mui * gi  * ri;
    pg  += pi  * gi  * ri;
    pp  += pi  * pi  * ri;
    pmu += pi  * mui * ri;
  }
  a = tau_scale + gg;
  b = mug - 2 * pg - eta * tau_scale;
  c = pp - pmu;
  return root_plus_stable_coeffs(a, b, c);
}

static const char *test_root_plus_equivalence(void) {
  /* Test with several hand-crafted input vectors */
  scs_int nm;
  scs_float tau_scale, eta, old_val, new_val, err;

  /* Test case 1: small vector */
  {
    scs_float g1[]  = {1.0, -2.0, 0.5};
    scs_float p1[]  = {0.3, 0.7, -0.1};
    scs_float mu1[] = {-0.5, 1.2, 0.8};
    scs_float r1[]  = {2.0, 0.5, 1.5};
    nm = 3;
    tau_scale = 1.0;
    eta = 0.5;
    old_val = root_plus_old(g1, p1, mu1, r1, nm, tau_scale, eta);
    new_val = root_plus_new(g1, p1, mu1, r1, nm, tau_scale, eta);
    err = ABS(old_val - new_val);
    scs_printf("root_plus case 1: old=%.12e new=%.12e err=%.2e\n",
               old_val, new_val, err);
    mu_assert("root_plus case 1 mismatch", err < 1e-12 * MAX(ABS(old_val), 1.));
  }

  /* Test case 2: larger vector with mixed signs */
  {
    scs_float g2[]  = {-0.1, 3.0, -2.5, 0.7, 1.1, -0.3, 0.9, -1.4};
    scs_float p2[]  = {0.5, -0.8, 1.2, -0.4, 0.6, 2.1, -1.0, 0.3};
    scs_float mu2[] = {1.0, -1.5, 0.3, 0.8, -0.2, 0.7, 1.3, -0.6};
    scs_float r2[]  = {0.1, 1.0, 3.0, 0.5, 2.0, 0.8, 1.5, 0.3};
    nm = 8;
    tau_scale = 2.5;
    eta = -0.3;
    old_val = root_plus_old(g2, p2, mu2, r2, nm, tau_scale, eta);
    new_val = root_plus_new(g2, p2, mu2, r2, nm, tau_scale, eta);
    err = ABS(old_val - new_val);
    scs_printf("root_plus case 2: old=%.12e new=%.12e err=%.2e\n",
               old_val, new_val, err);
    mu_assert("root_plus case 2 mismatch", err < 1e-12 * MAX(ABS(old_val), 1.));
  }

  /* Test case 3: large tau_scale dominates */
  {
    scs_float g3[]  = {0.01, -0.02};
    scs_float p3[]  = {100.0, -50.0};
    scs_float mu3[] = {200.0, 300.0};
    scs_float r3[]  = {1.0, 1.0};
    nm = 2;
    tau_scale = 1e6;
    eta = 1.0;
    old_val = root_plus_old(g3, p3, mu3, r3, nm, tau_scale, eta);
    new_val = root_plus_new(g3, p3, mu3, r3, nm, tau_scale, eta);
    err = ABS(old_val - new_val);
    scs_printf("root_plus case 3: old=%.12e new=%.12e err=%.2e\n",
               old_val, new_val, err);
    mu_assert("root_plus case 3 mismatch", err < 1e-10 * MAX(ABS(old_val), 1.));
  }

  /* Test case 4: near-zero discriminant */
  {
    scs_float g4[]  = {1.0, 0.0, 0.0, 0.0, 0.0};
    scs_float p4[]  = {0.0, 0.0, 0.0, 0.0, 0.0};
    scs_float mu4[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    scs_float r4[]  = {1.0, 1.0, 1.0, 1.0, 1.0};
    nm = 5;
    tau_scale = 1.0;
    eta = 0.0;
    old_val = root_plus_old(g4, p4, mu4, r4, nm, tau_scale, eta);
    new_val = root_plus_new(g4, p4, mu4, r4, nm, tau_scale, eta);
    err = ABS(old_val - new_val);
    scs_printf("root_plus case 4: old=%.12e new=%.12e err=%.2e\n",
               old_val, new_val, err);
    mu_assert("root_plus case 4 mismatch", err < 1e-14);
  }

  /* Test case 5: varying r weights spanning orders of magnitude */
  {
    scs_float g5[]  = {0.5, -1.3, 2.1, -0.7, 0.9, 1.1};
    scs_float p5[]  = {-0.2, 0.8, -1.5, 0.4, -0.6, 1.0};
    scs_float mu5[] = {0.3, -0.9, 0.6, 1.2, -0.8, 0.1};
    scs_float r5[]  = {1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6};
    nm = 6;
    tau_scale = 0.01;
    eta = 2.0;
    old_val = root_plus_old(g5, p5, mu5, r5, nm, tau_scale, eta);
    new_val = root_plus_new(g5, p5, mu5, r5, nm, tau_scale, eta);
    err = ABS(old_val - new_val);
    scs_printf("root_plus case 5: old=%.12e new=%.12e err=%.2e\n",
               old_val, new_val, err);
    mu_assert("root_plus case 5 mismatch", err < 1e-10 * MAX(ABS(old_val), 1.));
  }

  /* Test case 6: cancellation-prone direct formula */
  {
    scs_float stable_val, expected;
    stable_val = root_plus_stable_coeffs(1., 1e8, 1.);
    expected = -1e-8;
    err = ABS(stable_val - expected);
    scs_printf("root_plus case 6: stable=%.12e expected=%.12e err=%.2e\n",
               stable_val, expected, err);
    mu_assert("root_plus case 6 cancellation", err < 1e-14);
  }

  return 0;
}
