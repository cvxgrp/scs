#include "cones.h"
#include "linalg.h"
#include "scs_blas.h"
#include "scs_types.h"
#include "util.h"
#include "util_spectral_cones.h"

#include <stdlib.h>
#include <string.h>

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements code for projecting onto the ell1-norm cone.
 *
 * Last modified: 25 August 2024.
 */

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(axpy)(blas_int *n, const scs_float *a, const scs_float *x,
                blas_int *incx, scs_float *y, blas_int *incy);

scs_float BLAS(dot)(const blas_int *n, const scs_float *x, const blas_int *incx,
                    const scs_float *y, const blas_int *incy);

#ifdef __cplusplus
}
#endif

#ifdef SPECTRAL_DEBUG
static void compute_cone_residuals_ell1(const scs_float *tx, scs_float t0,
                                        const scs_float *x0, scs_int n,
                                        scs_float residuals[3]) {
  scs_float dual_res, pri_res, complementarity;
  scs_float dualt = tx[0] - t0;
  scs_float *dualx = malloc(n * sizeof(scs_float));
  blas_int int_n = n;
  scs_float negOne = -1.0;
  blas_int one = 1;
  scs_float ell1_norm = 0;
  scs_float inf_norm = 0;
  scs_float abs_val;
  const scs_float *xi;
  scs_int i;

  /* -------------------------------------
   *      Compute Lagrange multiplier.
   * (This function is not used in production so it is fine to allocate
   * memory here)
   * ------------------------------------- */
  memcpy(dualx, tx + 1, n * sizeof(scs_float));
  BLAS(axpy)(&int_n, &negOne, x0, &one, dualx, &one);

  /* ---------------------------------------
   *     Compute complementarity measure
   * --------------------------------------- */
  complementarity =
      tx[0] * dualt + BLAS(dot)(&int_n, dualx, &one, tx + 1, &one);

  /* -----------------------------------------------
   *   Compute primal feasibility measure
   * ----------------------------------------------- */
  for (xi = tx + 1; xi < tx + 1 + n; ++xi) {
    ell1_norm += fabs(*xi);
  }
  pri_res = ell1_norm - tx[0];

  /* ---------------------------------------
   *   Compute dual feasibility measure
   * --------------------------------------- */
  for (i = 0; i < n; ++i) {
    abs_val = fabs(dualx[i]);
    if (abs_val > inf_norm) {
      inf_norm = abs_val;
    }
  }
  dual_res = inf_norm - dualt;

  /* ------------------------------------------
   *  Assign result and free allocated memory
   * ------------------------------------------ */
  residuals[0] = dual_res;
  residuals[1] = pri_res;
  residuals[2] = complementarity;
  free(dualx);
}
#endif

/* Assumes that all components of x0 are nonnegative and
 * x0[0] >= x0[1] >= ... x0[n-1].
 *
 * Projects (t0, x0) onto {(t, x) : ||x||_1 <= t}. For sorted nonnegative
 * input the projection is proj_x[i] = max(x0[i] - theta, 0) and
 * proj_t = t0 + theta, where theta = (sum_{i<k} x0[i] - t0) / (k + 1) and
 * k is the largest index with x0[k-1] - theta(k) > 0 (Duchi et al. 2008,
 * "Efficient projections onto the l1-ball"). Scanning for the largest such
 * k is robust to repeated entries, where two-sided bracket searches can
 * fail to fire and produce points far outside the cone. */
scs_int ell1_cone_proj_sorted(scs_float t0, const scs_float *x0,
                              scs_float *proj, scs_int n) {
  scs_float xSum, theta, th;
  scs_int k, kk, i;

  if (-t0 >= x0[0]) {
    /* (t0, x0) lies in the polar cone: projection is the origin */
    memset(proj, 0, (n + 1) * sizeof(*x0));
    return 0;
  }

  xSum = 0;
  for (i = 0; i < n; ++i) {
    xSum += x0[i];
  }
  if (xSum <= t0) {
    /* (t0, x0) already lies in the cone: projection is the identity */
    proj[0] = t0;
    memcpy(proj + 1, x0, n * sizeof(*x0));
    return 0;
  }

  xSum = 0;
  theta = 0;
  k = 0;
  for (kk = 1; kk <= n; ++kk) {
    xSum += x0[kk - 1];
    th = (xSum - t0) / (kk + 1);
    if (x0[kk - 1] - th > 0) {
      k = kk;
      theta = th;
    } else {
      /* x0 is nonincreasing and theta(kk) nondecreasing: no later kk works */
      break;
    }
  }

  if (k == 0) {
    /* only reachable through rounding on polar-boundary inputs */
    memset(proj, 0, (n + 1) * sizeof(*x0));
    return 0;
  }

  proj[0] = MAX(t0 + theta, 0.);
  for (i = 0; i < k; ++i) {
    proj[i + 1] = x0[i] - theta;
  }
  memset(proj + 1 + k, 0, (n - k) * sizeof(*x0));
  return 0;
}

static void in_place_shuffle_ell1(scs_float *x, Value_index *work, scs_int n) {
  scs_int i;
  scs_int target_idx;
  scs_int temp_idx;
  scs_float temp_x;
  for (i = 0; i < n; ++i) {
    while (work[i].index != i) {
      /* Swap elements in `x` */
      target_idx = work[i].index;
      temp_x = x[i];
      x[i] = x[target_idx];
      x[target_idx] = temp_x;

      /* Swap indices in `idxs` to reflect the change */
      temp_idx = work[i].index;
      work[i].index = work[target_idx].index;
      work[target_idx].index = temp_idx;
    }
  }
}

int custom_cmp(const void *a, const void *b) {
  Value_index *elemA = (Value_index *)a;
  Value_index *elemB = (Value_index *)b;
  scs_float diff = fabs(elemB->value) - fabs(elemA->value);
  if (diff > 0) return 1;
  if (diff < 0) return -1;
  return 0;
}

void SCS(proj_ell_one)(scs_float *tx, scs_int n, ScsConeWork *c) {
  scs_float t0 = tx[0];
  scs_float *x0 = tx + 1;
  scs_float *proj = c->work_ell1_proj;
  Value_index *work = c->work_ell1;
  scs_int i;

  /* Preprocess vector so it is nonnegative and sorted */
  for (i = 0; i < n; ++i) {
    work[i].value = fabs(x0[i]);
    work[i].index = i;
  }

  qsort(work, n, sizeof(Value_index), custom_cmp);

  for (i = 0; i < n; ++i) {
    proj[i + 1] = work[i].value;
  }

  /* project preprocessed vector */
  ell1_cone_proj_sorted(t0, proj + 1, proj, n);

  /* recover original vector */
  in_place_shuffle_ell1(proj + 1, work, n);
  for (i = 0; i < n; i++) {
    proj[i + 1] = proj[i + 1] * (x0[i] >= 0 ? 1 : -1);
  }

  memcpy(tx, proj, (n + 1) * sizeof(*proj));
}
