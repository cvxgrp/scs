#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h"
#include "util_spectral_cones.h"

#define DUAL_FEAS_TOL 1e-2
#define PRI_FEAS_TOL 1e-2
#define COMP_TOL 1e-2
#define DUAL_T_THRESHOLD 1e-8
#define DUAL_X_THRESHOLD 1e-8

#define NEWTON_SUCCESS 1
#define IPM_VARIANT_0_SUCCESS 2
#define IPM_VARIANT_1_SUCCESS 3

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements a wrapper to the code that projects onto the logarithmic
 * cone.
 *
 * Last modified: 25 August 2024.
 */

// forward declare from log_cone_Newton.c
scs_int log_cone_Newton(scs_float t0, scs_float v0, const scs_float *x0,
                        scs_float *proj, scs_int n, scs_float *workspace,
                        Newton_stats *stats, bool *warm_start);

// forward declare form log_cone_IPM.c
scs_int log_cone_IPM(scs_float t0, scs_float v0, scs_float *x0, scs_float *u1,
                     scs_int n, scs_float *workspace, Newton_stats *stats,
                     scs_int variant);

// forward declare from this file
static void check_opt_cond_log_cone(const scs_float *tvx, scs_float t0,
                                    scs_float v0, const scs_float *x0,
                                    scs_int n, scs_float residuals[3],
                                    scs_float *dualx);

scs_int log_cone_proj_wrapper(scs_float t0, scs_float v0, scs_float *x0,
                              scs_float *proj, scs_int n, scs_float *workspace,
                              Newton_stats *stats, bool *warm_start) {
  scs_int status;
  // -----------------------------------------------------------------------
  // 1. Solve problem with Newton's method. From Lagrange multiplier theory
  //    we know that the t-component of the projection must always be greater
  //    than t0, so if proj[0] < t0 it means that Newton's method has converged
  //    to wrong point. We return if proj[0] >= t0 and the residuals suggest
  //    that we converged to the correct point.
  // 2. We warmstart Newton using the solution of the previous iteration
  //    except for the first iteration.
  // 3. The current implementation of the warmstart assumes that there is only
  //    one spectral matrix cone so 'proj' isn't overwritten.
  // ------------------------------------------------------------------------
  status = log_cone_Newton(t0, v0, x0, proj, n, workspace, stats, warm_start);

  if (proj[0] >= t0 - 0.1 * fabs(t0)) {
    check_opt_cond_log_cone(proj, t0, v0, x0, n, stats->residuals, workspace);
  } else {
    status = -1;
  }

  if (status == 0 && stats->residuals[0] < DUAL_FEAS_TOL &&
      stats->residuals[1] < PRI_FEAS_TOL &&
      fabs(stats->residuals[2]) < COMP_TOL) {
    stats->newton_success = NEWTON_SUCCESS;
    return 0;
  }

  // ------------------------------------------------------------------------
  // Solve problem with primal-dual IPM incorporating Mehrotra's correction
  // etc.
  // ------------------------------------------------------------------------
  status = log_cone_IPM(t0, v0, x0, proj, n, workspace, stats, 0);
  check_opt_cond_log_cone(proj, t0, v0, x0, n, stats->residuals, workspace);

  *warm_start = true; // next iteration Newton should be warmstarted

  if (status == 0 && stats->residuals[0] < DUAL_FEAS_TOL &&
      stats->residuals[1] < PRI_FEAS_TOL &&
      fabs(stats->residuals[2]) < COMP_TOL) {
    stats->newton_success = IPM_VARIANT_0_SUCCESS;
    return 0;
  }

  // ------------------------------------------------------------------------
  // Solve problem with primal-dual IPM without Mehrotra's correction
  // etc. (In all experiments in the paper by Cederberg and Boyd,
  // the IPM above solves the problem correctly so the code below never
  // runs. However, during development I ran into a (possibly pathological
  // case) where the first IPM fails. If this happens we run below another
  // version. This version of the IPM solves the problem that the first version
  // failed on.)
  // ------------------------------------------------------------------------
  status = log_cone_IPM(t0, v0, x0, proj, n, workspace, stats, 1);
  check_opt_cond_log_cone(proj, t0, v0, x0, n, stats->residuals, workspace);

  if (status == 0 && stats->residuals[0] < DUAL_FEAS_TOL &&
      stats->residuals[1] < PRI_FEAS_TOL &&
      fabs(stats->residuals[2]) < COMP_TOL) {
    stats->newton_success = IPM_VARIANT_1_SUCCESS;
    return 0;
  }

#ifdef SPECTRAL_DEBUG
  scs_printf("FAILURE: logarithmic cone projection");
  scs_printf("dual_res / pri_res / comp / iter: %.3e, %.3e, %.3e, %d\n",
             stats->residuals[0], stats->residuals[1], stats->residuals[2],
             stats->iter);

  scs_printf("Projecting the following point:\n %.10e, %.10e", t0, v0);
  for (scs_int i = 0; i < n; i++) {
    scs_printf(" %.10e", x0[i]);
  }
  scs_printf("\n");
#endif

  return -1;
}

// tvx = [t, v, x].
static void check_opt_cond_log_cone(const scs_float *tvx, scs_float t0,
                                    scs_float v0, const scs_float *x0,
                                    scs_int n, scs_float residuals[3],
                                    scs_float *dualx) {
  scs_float pri_res, dual_res, complementarity;

  // -------------------------------------------------------
  //           Compute Lagrange multiplier
  // -------------------------------------------------------
  scs_float dualt = tvx[0] - t0;
  if (fabs(dualt) < DUAL_T_THRESHOLD) {
    dualt = DUAL_T_THRESHOLD;
  }
  scs_float dualv = tvx[1] - v0;
  for (scs_int i = 0; i < n; ++i) {
    dualx[i] = tvx[i + 2] - x0[i];
    if (fabs(dualx[i]) < DUAL_X_THRESHOLD) {
      dualx[i] = DUAL_X_THRESHOLD;
    }
  }

  // ---------------------------------------------------------------------
  //          Compute complementarity measure
  // ---------------------------------------------------------------------
  complementarity =
      tvx[0] * dualt + tvx[1] * dualv + SCS(dot)(dualx, tvx + 2, n);

  // ----------------------------------------------------------------------
  //          Compute primal feasibility measure
  // ----------------------------------------------------------------------
  if (tvx[1] > 0 && is_pos(tvx + 2, n)) {
    pri_res = -tvx[1] * (sum_log(tvx + 2, n) - n * log(tvx[1])) - tvx[0];
  } else {
    pri_res = tvx[1] * tvx[1];

    if (tvx[0] < 0) {
      pri_res += tvx[0] * tvx[0];
    }

    for (const scs_float *xi = tvx + 2; xi < tvx + 2 + n; ++xi) {
      if (*xi < 0) {
        pri_res += (*xi) * (*xi);
      }
    }
  }

  // ----------------------------------------------------------------------
  //          Compute dual feasibility measure
  // ----------------------------------------------------------------------
  if (dualt > 0 && is_pos(dualx, n)) {
    dual_res = dualt * (n * log(dualt) - n - sum_log(dualx, n)) - dualv;
  } else {
    dual_res = dualt * dualt;

    if (dualv < 0) {
      dual_res += dualv * dualv;
    }

    for (scs_int i = 0; i < n; i++) {
      if (dualx[i] < 0) {
        dual_res += dualx[i] * dualx[i];
      }
    }
  }

  // --------------------------------------------------------
  //    Normalize the residuals and assign the result
  // --------------------------------------------------------
  scs_float dual_norm =
      sqrt(dualt * dualt + dualv * dualv + SCS(norm_sq)(dualx, n));
  scs_float pri_norm = SCS(norm_2)(tvx, n + 2);
  residuals[0] = dual_res / MAX(dual_norm, 1.0);
  residuals[1] = pri_res / MAX(pri_norm, 1.0);
  double scale = MAX(pri_norm, 1.0);
  residuals[2] = complementarity / MAX(scale, dual_norm);
}

