#include "cones.h"
// #include "scs.h"
#include "linalg.h"
#include "scs_blas.h"
#include "scs_types.h"
#include "util.h" // just for timer
#include "util_spectral_cones.h"
#include <stdlib.h> // qsort

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

  // -------------------------------------
  //      Compute Lagrange multiplier.
  // (This function is not used in production so it is fine to allocate
  // memory here)
  // -------------------------------------
  scs_float dualt = tx[0] - t0;
  scs_float *dualx = malloc(n * sizeof(scs_float));
  memcpy(dualx, tx + 1, n * sizeof(scs_float));
  blas_int int_n = n;
  scs_float negOne = -1.0;
  blas_int one = 1;
  BLAS(axpy)(&int_n, &negOne, x0, &one, dualx, &one);

  // ---------------------------------------
  //     Compute complementarity measure
  // ---------------------------------------
  complementarity =
      tx[0] * dualt + BLAS(dot)(&int_n, dualx, &one, tx + 1, &one);

  // -----------------------------------------------
  //   Compute primal feasibility measure
  // -----------------------------------------------
  scs_float ell1_norm = 0;
  for (const scs_float *xi = tx + 1; xi < tx + 1 + n; ++xi) {
    ell1_norm += fabs(*xi);
  }
  pri_res = ell1_norm - tx[0];

  // ---------------------------------------
  //   Compute dual feasibility measure
  // ---------------------------------------
  scs_float inf_norm = 0;
  for (scs_int i = 0; i < n; ++i) {
    scs_float abs_val = fabs(dualx[i]);
    if (abs_val > inf_norm) {
      inf_norm = abs_val;
    }
  }
  dual_res = inf_norm - dualt;

  // ------------------------------------------
  //  Assign result and free allocated memory
  // ------------------------------------------
  residuals[0] = dual_res;
  residuals[1] = pri_res;
  residuals[2] = complementarity;
  free(dualx);
}
#endif

// Asssumes that all components of x0 are positive and
// x0[0] >= x0[1] >= ... x0[n-1].
scs_int ell1_cone_proj_sorted(scs_float t0, const scs_float *x0,
                              scs_float *proj, scs_int n) {
  if (-t0 >= x0[0]) {
    memset(proj, 0, (n + 1) * sizeof(*x0));
    return 0;
  }

  // -------------------------------------------
  //            Find the value on k
  // -------------------------------------------

  // check if k = 0 suffices
  if (-t0 >= x0[0]) {
    memset(proj, 0, (n + 1) * sizeof(*x0));
    return 0;
  }

  scs_float xSum = 0;
  scs_float tempSum = 0;
  int k = -1;
  for (scs_int kk = 1; kk < n; ++kk) {
    xSum += x0[kk - 1];
    tempSum = (-t0 + xSum) / (kk + 1);

    if (x0[kk - 1] > tempSum && x0[kk] <= tempSum) {
      k = (int)kk;
      break;
    }
  }

  if (k == -1) {
    k = n;
    xSum += x0[n - 1];
  }

  // ---------------------------------------------
  //                Execute projection
  // ---------------------------------------------
  proj[0] = -t0 + xSum;

  if (proj[0] > 0) {
    proj[0] = t0 + proj[0] / (k + 1);
  } else {
    proj[0] = t0;
  }

  memcpy(proj + 1, x0, k * sizeof(*x0));
  scs_float diff = proj[0] - t0;
  for (int i = 1; i < k + 1; i++) {
    proj[i] -= diff;
  }
  memset(proj + 1 + k, 0, (n - k) * sizeof(*x0));

#ifdef SPECTRAL_DEBUG
  //-------------------------------------------------------------------------
  //               Check residuals - not needed in production
  //-------------------------------------------------------------------------
  scs_float residuals[3];
  compute_cone_residuals_ell1(proj, t0, x0, n, residuals);

  if (residuals[0] > 1e-8 || residuals[1] > 1e-8 || residuals[2] > 1e-8) {
    scs_printf("WARN: something is wrong in nuclear norm cone projection.\n");
    scs_printf("dual_res / primal_res / comp : %.3e, %.3e, %.3e\n",
               residuals[0], residuals[1], residuals[2]);
    return -1;
  }
#endif

  return 0;
}

static void in_place_shuffle_ell1(scs_float *x, Value_index *work, scs_int n) {
  for (scs_int i = 0; i < n; ++i) {
    while (work[i].index != i) {
      // Swap elements in `x`
      scs_int target_idx = work[i].index;
      scs_float temp_x = x[i];
      x[i] = x[target_idx];
      x[target_idx] = temp_x;

      // Swap indices in `idxs` to reflect the change
      scs_int temp_idx = work[i].index;
      work[i].index = work[target_idx].index;
      work[target_idx].index = temp_idx;
    }
  }
}

int custom_cmp(const void *a, const void *b) {
  Value_index *elemA = (Value_index *)a;
  Value_index *elemB = (Value_index *)b;
  return fabs(elemB->value) - fabs(elemA->value) > 0 ? 1 : -1;
}

void SCS(proj_ell_one)(scs_float *tx, scs_int n, ScsConeWork *c) {
  scs_float t0 = tx[0];
  scs_float *x0 = tx + 1;
  scs_float *proj = c->work_ell1_proj;
  Value_index *work = c->work_ell1;
  // ------------------------------------------------------
  //     Preprocess vector so it is nonnegative and sorted
  // ------------------------------------------------------
  for (scs_int i = 0; i < n; ++i) {
    work[i].value = fabs(x0[i]);
    work[i].index = i;
  }

  qsort(work, n, sizeof(Value_index), custom_cmp);

  for (scs_int i = 0; i < n; ++i) {
    proj[i + 1] = work[i].value;
  }

  // ------------------------------------------
  //        project preprocessed vector
  // ------------------------------------------
  ell1_cone_proj_sorted(t0, proj + 1, proj, n);

  // -------------------------------------------
  //        recover original vector
  // -------------------------------------------
  in_place_shuffle_ell1(proj + 1, work, n);
  for (scs_int i = 0; i < n; i++) {
    proj[i + 1] = proj[i + 1] * (x0[i] >= 0 ? 1 : -1);
  }

  memcpy(tx, proj, (n + 1) * sizeof(*proj));
}
