#include "cones.h"
// #include "scs.h"
#include "linalg.h"
#include "scs_blas.h"
#include "scs_types.h"
#include "util_spectral_cones.h"
#include <stdlib.h> // qsort

#define TOL_LARGEST_CONE 1e-9

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements code for projecting onto the sum-of-largest cone.
 * It assumes that the input is sorted. If you need code that does not make
 * this assumption, please reach out.
 *
 * Last modified: 25 August 2024.
 */

#ifdef SPECTRAL_DEBUG
static void compute_cone_residuals(scs_float t, const scs_float *x,
                                   scs_float t0, const scs_float *x0,
                                   scs_float residuals[3], scs_int n,
                                   scs_int k);
#endif

scs_int assert_sorted(scs_float *x, scs_int n) {
  for (scs_int i = 0; i < n - 1; ++i) {
    if (x[i] < x[i + 1]) {
      return -1;
    }
  }
  return 0;
}

scs_int proj_sum_largest_cone_sorted(scs_float *t, scs_float *x, scs_int n,
                                     scs_int k) {
#ifdef SPECTRAL_DEBUG
  scs_float t00 = *t;
  scs_float *x0 = scs_malloc(n * sizeof(*x0));
  memcpy(x0, x, n * sizeof(*x0));
  scs_int status = assert_sorted(x, n);
  if (status < 0) {
    scs_printf("NOT SORTED! \n");
    return status;
  }
#endif

  // -------------------------------
  //    Initialize state variables
  // -------------------------------
  assert(k < n && k > 0);
  scs_int nu = k, nt = 0;
  scs_float eta = 0.0;
  scs_float t0 = *t;
  scs_float S = x[0];
  for (scs_int i = 1; i < k; ++i) {
    S += x[i];
  }

  scs_float a_u = x[nu - 1], a_t = x[nu];

  // ---------------------------------
  //           main loop
  // ---------------------------------
  scs_float ratio, s1, s3, s;
  while (S > *t + TOL_LARGEST_CONE) {
    ratio = (nu == k) ? 1.0 : (scs_float)(nt) / (k - nu);

    // ------------------------------------------------------
    //                 compute step size
    // ------------------------------------------------------
    s1 = (nu == k) ? a_u - a_t : (a_u - a_t) / (ratio - 1);
    s3 = (S - *t) / (ratio * (nu + 1) + (k - nu));
    s = (nu == 0) ? s3 : MIN(s3, s1);

    if (!(nu + nt == n || nt == 0)) {
      scs_float val = a_t - x[nu + nt];
      s = MIN(s, val);
    }

    // --------------------------
    //        update state
    // --------------------------
    eta += s * ratio;
    S -= s * (ratio * nu + k - nu);
    *t = t0 + eta;

    if (nt > 0) {
      a_t -= s;
    }

    if (nu != 0 && s == s1) {
      nu -= 1;
    } else {
      assert(s != s1);
    }

    if (nu > 0) {
      a_u = x[nu - 1] - eta;
    }

    nt = (nt == 0) ? 2 : nt + 1;
  }

  nt -= 1;

  // update projection
  scs_int i;
  for (i = 0; i < nu; ++i) {
    x[i] -= eta;
  }

  for (i = nu; i < nu + nt; ++i) {
    x[i] = a_t;
  }

#ifdef SPECTRAL_DEBUG
  scs_float residuals[3];
  compute_cone_residuals(*t, x, t00, x0, residuals, n, k);
  scs_free(x0);
  if (residuals[0] > 1e-10 || residuals[1] > 1e-10 ||
      fabs(residuals[2]) > 1e-10) {
    return -1;
  }
#endif

  return 0;
}

scs_int cmp_desc(const void *a, const void *b) {
  scs_float da = *(const scs_float *)a;
  scs_float db = *(const scs_float *)b;
  if (da < db)
    return 1;
  if (da > db)
    return -1;
  return 0;
}

#ifdef SPECTRAL_DEBUG
// this function is not used in production so fine to allocate memory
static scs_float sum_largest_val(const scs_float *x, scs_int n, scs_int k) {
  scs_float *x_temp = scs_malloc(n * sizeof(*x));
  memcpy(x_temp, x, n * sizeof(*x));
  qsort(x_temp, n, sizeof(*x_temp), cmp_desc);

  scs_float val = 0.0;
  for (scs_int i = 0; i < k; ++i) {
    val += x_temp[i];
  }
  scs_free(x_temp);
  return val;
}

static void compute_cone_residuals(scs_float t, const scs_float *x,
                                   scs_float t0, const scs_float *x0,
                                   scs_float residuals[3], scs_int n,
                                   scs_int k) {
  scs_float lmbda_t = t - t0;
  scs_float *lmbda_x = scs_malloc(n * sizeof(*x));
  scs_float sum_lmbda_x = 0.0;
  for (scs_int i = 0; i < n; ++i) {
    lmbda_x[i] = x[i] - x0[i];
    sum_lmbda_x += lmbda_x[i];
  }

  scs_float comp = lmbda_t * t + SCS(dot)(lmbda_x, x, n);
  scs_float pri_res = sum_largest_val(x, n, k) - t;

  scs_float dual_res = fabs(sum_lmbda_x + lmbda_t * k);
  dual_res *= dual_res;

  for (scs_int i = 0; i < n; ++i) {
    if (lmbda_x[i] > 0) {
      dual_res += lmbda_x[i] * lmbda_x[i];
    }

    if (lmbda_x[i] + lmbda_t < 0) {
      dual_res += (lmbda_x[i] + lmbda_t) * (lmbda_x[i] + lmbda_t);
    }
  }

  residuals[0] = dual_res;
  residuals[1] = pri_res;
  residuals[2] = comp;

  scs_free(lmbda_x);
}
#endif

