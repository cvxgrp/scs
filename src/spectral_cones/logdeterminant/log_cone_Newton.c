#include "glbopts.h" // for scs_printf
#include "linalg.h"
#include "scs_blas.h"
#include "util_spectral_cones.h"
#include <string.h> // for memcpy

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements Newton's method for projecting onto the logarithmic
 * cone.
 *
 * Last modified: 25 August 2024.
 */

#define LINESEARCH_RELATIVE_TOL 1e-14
#define MIN_INIT_LOG_CONE 1
#define MIN_DENOMINATOR 1e-14
#define MIN_X 1e-17
#define MIN_FLOAT MIN_X / 2
#define MIN_V 1e-14

#define MAX_ITER_NEWTON 100
#define ALPHA_NEWTON 0.01
#define BETA_NEWTON 0.8
#define TOL_NEWTON 1e-12
#define MAX_GRAD_STEPS 5

#define TERMINATE_DUE_TO_ZEROS -5
#define MAX_GRAD_STEPS_REACHED -6

// the CALLER of this function must make sure that the arguments belong to the
// domain
static scs_float obj_val(const scs_float *u, scs_float t0, scs_float v0,
                         const scs_float *x0, scs_int n) {
  scs_float v = u[0];
  const scs_float *x = u + 1;

  assert(v > 0 && min_vec(x, n) > 0);

  scs_float sx = -(v * sum_log(x, n) - n * v * log(v));
  scs_float obj = 0.5 * (sx - t0) * (sx - t0) + 0.5 * (v - v0) * (v - v0);
  for (scs_int i = 0; i < n; ++i) {
    obj += 0.5 * (x[i] - x0[i]) * (x[i] - x0[i]);
  }
  return obj;
}

scs_int log_cone_Newton(scs_float t0, scs_float v0, const scs_float *x0,
                        scs_float *u, scs_int n, scs_float *workspace,
                        Newton_stats *stats, bool *warm_start) {
  scs_float *t = u;
  scs_float *v = u + 1;
  scs_float *x = u + 2;
  scs_int n_plus_one = n + 1;

  // -------------------------------------------------------------------------
  //     check cone membership and special cases with analytic projections
  // -------------------------------------------------------------------------
  bool ix_x0_non_neg = true;
  bool is_x0_pos = true;
  scs_int i = 0;
  while (ix_x0_non_neg || is_x0_pos) {
    if (x[i] < 0) {
      ix_x0_non_neg = false;
      is_x0_pos = false;
    } else if (x[i] == 0) {
      is_x0_pos = false;
    }

    if (i == n - 1) {
      break;
    }
    ++i;
  }

  // if (t0, v0, x0) belongs to cone
  if ((v0 > 0 && is_x0_pos && -v0 * (sum_log(x0, n) - n * log(v0)) <= t0) ||
      (v0 == 0 && ix_x0_non_neg && t0 >= 0)) {
    *t = t0;
    *v = v0;
    memcpy(x, x0, sizeof(*x0) * n);
    stats->iter = IN_CONE;
    *warm_start = true;
    return 0;
  }

  // if (t0, v0, x0) belongs to negative dual cone
  if (t0 < 0 && is_negative(x0, n)) {
    scs_float sum = -n;
    for (i = 0; i < n; ++i) {
      sum -= log(x0[i] / t0);
    }
    sum *= t0;

    if (v0 <= sum) {
      memset(u, 0, (n + 2) * sizeof(*x0));
      stats->iter = IN_NEGATIVE_DUAL_CONE;
      // if 0 is the solution we should not use it to warmstart the next
      // iteration
      *warm_start = false;
      return 0;
    }
  }

  // special case with analytic solution
  if (v0 <= 0 && t0 >= 0) {
    *t = t0;
    *v = 0;
    non_neg_proj(x0, x, n);
    stats->iter = ANALYTICAL_SOL;
    *warm_start = true;
    return 0;
  }

  // ----------------------------------------------------------------------
  // if 'warm_start' is false we initialize in the point
  // (v, x) = (max(v0, MIN_INIT_LOG_CONE), max(x0, MIN_INIT_LOG_CONE)),
  // otherwise it is assumed that 'proj' has already been
  // initialized / warmstarted.
  // ----------------------------------------------------------------------
  if (!(*warm_start)) {
    *v = (v0 > MIN_INIT_LOG_CONE) ? v0 : MIN_INIT_LOG_CONE;

    for (scs_int i = 0; i < n; ++i) {
      x[i] = (x0[i] > MIN_INIT_LOG_CONE) ? x0[i] : MIN_INIT_LOG_CONE;
    }
  }

  scs_float obj_old = obj_val(u + 1, t0, v0, x0, n);

  // -----------------------------
  //      parse workspace
  // -----------------------------
  scs_float *grad = workspace;
  scs_float *d = grad + n_plus_one;
  scs_float *w = d + n_plus_one;
  scs_float *du = w + n_plus_one;
  scs_float *temp1 = du + n_plus_one;
  scs_float *u_new = du + n_plus_one + 1;

  size_t num_grad_steps = 0;
  scs_int iter = 0;
  scs_float newton_decrement = 0;
  for (iter = 1; iter <= MAX_ITER_NEWTON; ++iter) {
    // A small value on v indicates that Newton's method converges to the
    // origin. In this case we should abort and apply an IPM.
    if (*v < MIN_V) {
      stats->iter = TERMINATE_DUE_TO_ZEROS;
      return -1;
    }

    // -------------------------------------------------------------------
    // To avoid pathological cases where some components of
    // x approach 0 we use a minimum threshold.
    // -------------------------------------------------------------------
    for (scs_int i = 0; i < n; i++) {
      x[i] = MAX(x[i], MIN_X);
    }

    // ----------------------------------------------------------------
    //                  compute gradient and Hessian
    // ----------------------------------------------------------------
    assert(*v > MIN_FLOAT && min_vec(x, n) > MIN_FLOAT);
    scs_float temp0 = -sum_log(x, n) + n * log(*v);
    scs_float a = (*v) * temp0 - t0;
    scs_float c = temp0 + n;

    grad[0] = a * c + (*v) - v0;
    scs_float v_inv = 1 / (*v);
    d[0] = 1 + a * (-a * (v_inv * v_inv) + n * v_inv - 2 * c * v_inv);
    w[0] = -(a + (*v) * c) * v_inv;
    scs_float av = a * (*v);

    for (i = 1; i < n + 1; ++i) {
      assert(x[i - 1] > 0);
      scs_float x_inv = 1 / x[i - 1];
      grad[i] = -av * x_inv + x[i - 1] - x0[i - 1];
      d[i] = 1 + av * (x_inv * x_inv);
      w[i] = (*v) * x_inv;
    }

    // ----------------------------------------------------------------------
    //    Solve for Newton step. I have seen a scenario when the variable
    //   'denominator' becomes 0. This occurred when Newton's method
    //    converged to the origin and the IPM was necessary
    //    to converge to the correct point. We therefore abort when the
    //    denominator becomes very close to 0, since this may indicate
    //    that Newton's method is converging to the wrong point and in this
    //    case it is necessary to apply an IPM.
    // ----------------------------------------------------------------------
    scs_float nominator = 0;
    scs_float denominator = 1;
    for (i = 0; i < n + 1; ++i) {
      assert(fabs(d[i]) > MIN_FLOAT);
      du[i] = -grad[i] / d[i];
      nominator += w[i] * du[i];
      temp1[i] = w[i] / d[i];
      denominator += w[i] * temp1[i];
    }

    if (fabs(denominator) < MIN_DENOMINATOR) {
      stats->iter = TERMINATE_DUE_TO_ZEROS;
      return -1;
    }

    scs_float ratio = -nominator / denominator;
    SCS(add_scaled_array)(du, temp1, n_plus_one, ratio);

    // --------------------------------------------------------------------
    // if the Newton direction is not descent we use the negative gradient
    // as the search direction, provided that it hasn't been used many times
    // before.
    // --------------------------------------------------------------------
    scs_float dirDer = SCS(dot)(grad, du, n_plus_one);
    if (dirDer > 0) {
      if (num_grad_steps >= MAX_GRAD_STEPS) {
        stats->iter = MAX_GRAD_STEPS_REACHED;
        return -1;
      }

      num_grad_steps += 1;
      dirDer = 0;

      for (i = 0; i < n + 1; ++i) {
        du[i] = -grad[i];
        dirDer -= grad[i] * grad[i];
      }
    }

    // --------------------------------------------------------
    // compute Newton decrement and check termination criteria
    // -------------------------------------------------------
    newton_decrement = -dirDer;
    if (newton_decrement <= 2 * TOL_NEWTON) {
      break;
    }

    // --------------------------------------------------
    // find largest step size with respect to the domain
    // --------------------------------------------------
    scs_float step_size = 1;
    for (i = 0; i < n + 1; i++) {
      if (du[i] < 0) {
        scs_float max_step = -0.99 * u[i + 1] / du[i];
        step_size = MIN(step_size, max_step);
      }
    }

    // -------------------------------------------------
    // backtracking line search. First two lines do
    // u_new = u + t * du;
    // -------------------------------------------------
    memcpy(u_new + 1, u + 1, n_plus_one * sizeof(*u));
    SCS(add_scaled_array)(u_new + 1, du, n_plus_one, step_size);
    scs_float new_obj = obj_val(u_new + 1, t0, v0, x0, n);
    while ((1 - LINESEARCH_RELATIVE_TOL) * new_obj >
           obj_old + ALPHA_NEWTON * step_size * dirDer) {
      step_size *= BETA_NEWTON;
      memcpy(u_new + 1, u + 1, n_plus_one * sizeof(*u));
      SCS(add_scaled_array)(u_new + 1, du, n_plus_one, step_size);
      new_obj = obj_val(u_new + 1, t0, v0, x0, n);
    }

    obj_old = new_obj;
    memcpy(u + 1, u_new + 1, n_plus_one * sizeof(*u));
  }

  assert(min_vec(x, n) > MIN_FLOAT && *v > MIN_FLOAT);
  *t = -(*v) * (sum_log(x, n) - n * log(*v));

  *warm_start = true;
  stats->iter = iter;
  return 0;
}
