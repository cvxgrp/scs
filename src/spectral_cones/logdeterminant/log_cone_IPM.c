#include "glbopts.h" // for scs_printf
#include "linalg.h"
#include "scs_blas.h"
#include "util_spectral_cones.h" // for Newton structs
#include <string.h>              // for memcpy

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements a primal-dual IPM for projecting onto the logarithmic
 * cone.
 *
 * Last modified: 25 August 2024.
 */

#define FEASTOL_IPM 1e-7
#define ABSTOL_IPM 1e-7
#define RELTOL_IPM 1e-6
#define MAX_ITER_IPM 100
#define BETA_IPM 0.5
#define STEP_IPM 0.99
#define ALPHA_IPM 0.01
#define MAX_RELAXED_ITERS 8

blas_int bi_one = 1;
scs_float d_one = 1.0;
scs_float d_minus_one = -1.0;
blas_int bi_three = 3;

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(gemv)(const char *trans, const blas_int *m, const blas_int *n,
                const scs_float *alpha, const scs_float *a, const blas_int *lda,
                const scs_float *x, const blas_int *incx, const scs_float *beta,
                scs_float *y, const blas_int *incy);

#ifdef __cplusplus
}
#endif

/* Evaluates the oracle f(u) = (f0(u), f1(u), f2(u)) where u = [t, v, x, r]
 *  and
 *  f0 = 0.5 * (t - t0)^2 + 0.5 * (v - v0)^2 + 0.5 * ||x - x0||^2  - r
 *  f1 = - v * sum(log(x/v)) - t
 *  f2 = -v
 *
 *  Note that u1 = [t, v, x].
 */
static void f_oracle(const scs_float *u1, scs_float r, scs_float t0,
                     scs_float v0, const scs_float *x0, const scs_float *x_inv,
                     scs_int n, scs_float *f_u, scs_float *grad_f0,
                     scs_float *grad_f1) {
  // compute grad_f0[0, 1] = t - t0, v - v0, grad_f1[2:-1] = x - x0
  grad_f0[0] = u1[0] - t0;
  grad_f0[1] = u1[1] - v0;
  memcpy(grad_f0 + 2, u1 + 2, n * sizeof(*grad_f0));
  SCS(add_scaled_array)(grad_f0 + 2, x0, n, -1.0);
  grad_f0[n + 2] = -1.0;

  // compute f_u
  scs_float sum_log_xv = sum_log(u1 + 2, n) - n * log(u1[1]);
  scs_float norm_2 = SCS(norm_2)(grad_f0 + 2, n);
  f_u[0] = 0.5 * (grad_f0[0] * grad_f0[0] + grad_f0[1] * grad_f0[1] +
                  norm_2 * norm_2) -
           r;
  f_u[1] = -u1[1] * sum_log_xv - u1[0];
  f_u[2] = -u1[1];

  // compute grad_f1[0] = -1, grad_f1[1] = n - sum(log(x/v)),
  // grad[2:-1] = - v / x, grad[-1] = 0
  grad_f1[0] = -1.0;
  grad_f1[1] = n - sum_log_xv;
  memcpy(grad_f1 + 2, x_inv, n * sizeof(*x_inv));
  SCS(scale_array)(grad_f1 + 2, -u1[1], n);
  grad_f1[n + 2] = 0.0;
}

// u1 = [t, v, x]
static scs_float find_max_step_size(const scs_float *u1, const scs_float *z,
                                    const scs_float *s, const scs_float *du1,
                                    const scs_float *dz, const scs_float *ds,
                                    scs_int n) {
  scs_float step_max = 10.0;
  scs_int i;

  // find maximum step size with respect to z and s
  for (i = 0; i < 3; ++i) {
    if (dz[i] < 0) {
      step_max = MIN(step_max, -z[i] / dz[i]);
    }

    if (ds[i] < 0) {
      step_max = MIN(step_max, -s[i] / ds[i]);
    }
  }

  scs_float step_max_dom = 10;

  // find maximum step size with respect to the domain of f
  for (i = 1; i < n + 2; i++) {
    if (du1[i] < 0) {
      step_max_dom = MIN(step_max_dom, -u1[i] / du1[i]);
    }
  }

  scs_float step_size = MIN(1.0, STEP_IPM * step_max);

  while (step_size > step_max_dom) {
    step_size *= BETA_IPM;
  }

  return step_size;
}

typedef struct {
  scs_float *temp1;
  scs_float *temp2;
  scs_float *GinvC;
  scs_float *R;
  scs_float *g0g1;
  scs_float coeff;
  scs_float *GinvRes;
  scs_float *CTGinvRes;
  scs_float *RinvCTGinvRes;
  scs_float *residual;
  scs_float *rhs_reduced;
  scs_float *bnew;
  scs_float *CCTdu;

} KKT_solver_workspace;

// Precomputations before solving the KKT system. You can think of this
// as factoring the KKT matrix.
static void KKT_precompute(scs_float v, const scs_float *x_inv,
                           const scs_float *z, const scs_float *w, scs_int n,
                           KKT_solver_workspace *KKT_work) {
  scs_int u_dim = n + 3;
  scs_float *temp1 = KKT_work->temp1;
  scs_float *temp2 = KKT_work->temp2;
  scs_float *GinvC = KKT_work->GinvC;
  scs_float *g0g1 = KKT_work->g0g1;
  scs_float *R = KKT_work->R;

  scs_float a = z[0] + 1 / (w[2] * w[2]) + z[1] * n / v;
  scs_float coeff = 0.0;
  for (scs_int i = 0; i < n; i++) {
    temp1[i] = z[0] + (z[1] * v) * (x_inv[i] * x_inv[i]);
    temp2[i] = x_inv[i] / temp1[i];
    coeff += (x_inv[i] * x_inv[i]) / temp1[i];
  }
  coeff = a - z[1] * z[1] * coeff;
  KKT_work->coeff = coeff;

  // --------------------------------------------------------------------
  //                      compute G \ g0
  // --------------------------------------------------------------------
  GinvC[0] = g0g1[0] / z[0];
  GinvC[1] = (g0g1[1] + z[1] * SCS(dot)(g0g1 + 2, temp2, n)) / coeff;
  for (scs_int i = 2; i < n + 2; i++) {
    GinvC[i] = (g0g1[i] + (z[1] * GinvC[1]) * x_inv[i - 2]) / temp1[i - 2];
  }
  GinvC[n + 2] = -g0g1[n + 2];

  // --------------------------------------------------------------------
  //                      compute G \ g1
  // --------------------------------------------------------------------
  GinvC[n + 3] = g0g1[n + 3] / z[0];
  GinvC[n + 4] =
      (g0g1[n + 4] + z[1] * SCS(dot)(g0g1 + n + 5, temp2, n)) / coeff;
  for (scs_int i = n + 5; i < 2 * n + 5; i++) {
    GinvC[i] =
        (g0g1[i] + (z[1] * GinvC[n + 4]) * x_inv[i - n - 5]) / temp1[i - n - 5];
  }
  GinvC[2 * n + 5] = -g0g1[2 * n + 5];

  // compute G \ ([0, 0, ..., 0, 1])
  GinvC[3 * n + 8] = -1.0;

  // ------------------------------------------------------------------------
  //  compute R = I + C.T @ GinvC. This matrix is always reverse triangular.
  // ------------------------------------------------------------------------
  R[0] = 1 + SCS(dot)(g0g1, GinvC, u_dim);
  R[3] = SCS(dot)(g0g1, GinvC + u_dim, u_dim);
  R[6] = SCS(dot)(g0g1, GinvC + 2 * u_dim, u_dim);
  R[1] = SCS(dot)(g0g1 + u_dim, GinvC, u_dim);
  R[4] = 1 + SCS(dot)(g0g1 + u_dim, GinvC + u_dim, u_dim);
  R[2] = GinvC[n + 2];
}

static void KKT_solve(const scs_float *z, const scs_float *w,
                      const scs_float *x_inv, const scs_float v,
                      const scs_float *lmbda, scs_int n, const scs_float *rhs1,
                      const scs_float *rhs2, KKT_solver_workspace *KKT_work,
                      scs_float *du1, scs_float *dr, scs_float *dz,
                      scs_float *ds) {
  blas_int u_dim = n + 3;
  scs_int u1_dim = n + 2;

  scs_float *temp1 = KKT_work->temp1;
  scs_float *temp2 = KKT_work->temp2;
  scs_float *GinvC = KKT_work->GinvC;
  scs_float *R = KKT_work->R;
  scs_float coeff = KKT_work->coeff;
  scs_float *GinvRes = KKT_work->GinvRes;
  scs_float *CTGinvRes = KKT_work->CTGinvRes;
  scs_float *RinvCTGinvRes = KKT_work->RinvCTGinvRes;
  scs_float *residual = KKT_work->residual;
  scs_float *rhs_reduced = KKT_work->rhs_reduced;
  scs_float *bnew = KKT_work->bnew;
  scs_float *g0g1 = KKT_work->g0g1;
  scs_float *CCTdu = KKT_work->CCTdu;

  // -------------------------------------------------------------------------
  //                      prepare RHS
  // -------------------------------------------------------------------------
  memcpy(rhs_reduced, rhs1, (n + 6) * sizeof(*rhs1));
  rhs_reduced[u_dim] -= w[0] * (rhs2[0] / lmbda[0]);
  rhs_reduced[u_dim + 1] -= w[1] * (rhs2[1] / lmbda[1]);
  rhs_reduced[u_dim + 2] -= w[2] * (rhs2[2] / lmbda[2]);
  memcpy(bnew, rhs_reduced, u_dim * sizeof(*rhs_reduced));
  scs_float scale = rhs_reduced[n + 3] / w[0];
  SCS(add_scaled_array)(bnew, g0g1, u_dim, scale);
  scale = rhs_reduced[n + 4] / w[1];
  SCS(add_scaled_array)(bnew, g0g1 + u_dim, u_dim, scale);
  bnew[1] -= (rhs_reduced[n + 5] / (w[2] * w[2]));

  // -------------------------------------------------------------------------
  //                 solve system and apply iterative refinement
  // -------------------------------------------------------------------------
  memcpy(residual, bnew, u_dim * sizeof(*bnew));
  memset(du1, 0, u1_dim * sizeof(*du1));
  *dr = 0;
  const scs_int num_ref = 3;
  for (scs_int i = 0; i < num_ref; ++i) {
    // --------------------------------
    // solve (G + C @ C.T) d = residual
    // --------------------------------
    GinvRes[0] = residual[0] / z[0];
    GinvRes[1] =
        (residual[1] + z[1] * SCS(dot)(residual + 2, temp2, n)) / coeff;
    for (scs_int i = 2; i < n + 2; ++i) {
      GinvRes[i] =
          (residual[i] + (z[1] * GinvRes[1]) * x_inv[i - 2]) / temp1[i - 2];
    }

    GinvRes[n + 2] = -residual[n + 2];
    CTGinvRes[0] = SCS(dot)(g0g1, GinvRes, u_dim);
    CTGinvRes[1] = SCS(dot)(g0g1 + u_dim, GinvRes, u_dim);
    CTGinvRes[2] = GinvRes[n + 2];
    RinvCTGinvRes[0] = CTGinvRes[2] / R[2];
    RinvCTGinvRes[1] = (CTGinvRes[1] - R[1] * RinvCTGinvRes[0]) / R[4];
    RinvCTGinvRes[2] =
        (CTGinvRes[0] - R[0] * RinvCTGinvRes[0] - R[3] * RinvCTGinvRes[1]) /
        R[6];

    // --------------------------------------------------------------------
    //  Here we want to apply the correction d: du = du + d, where
    //  d = GinvRes - GinvC @ RinvCTGinvRes. Note that we expect the
    //  correction d to be very small, so we *must* compute d first and
    //  then add it to du, instead of computing it as
    //  du = du + GinvRes, du = du - GinvC @ RinvCTGinvRes.
    //  In other words, the following code will result in bugs due to
    //  numerics:
    //  daxpy_(&u_dim, &d_one, GinvRes, &bi_one, du, &bi_one);
    //  char trans = 'N';
    //  dgemv_(&trans, &u_dim, &i_three, &d_minus_one, GinvC, &u_dim,
    //         RinvCTGinvRes, &bi_one, &d_one, du, &bi_one);
    // --------------------------------------------------------------------

    char trans = 'N';
    BLAS(gemv)(&trans, &u_dim, &bi_three, &d_minus_one, GinvC, &u_dim,
               RinvCTGinvRes, &bi_one, &d_one, GinvRes, &bi_one);
    SCS(add_scaled_array)(du1, GinvRes, u1_dim, 1.0);
    *dr += GinvRes[n + 2];

    // -----------------------------------
    //        compute new residual
    // -----------------------------------
    scs_float *Gdu = GinvRes;
    if (i != num_ref - 1) {
      Gdu[0] = z[0] * du1[0];
      Gdu[1] = ((z[0] + 1 / (w[2] * w[2])) * du1[1] +
                z[1] * (n / v * du1[1] - SCS(dot)(x_inv, du1 + 2, n)));
      for (scs_int ii = 2; ii < n + 2; ++ii) {
        Gdu[ii] = (z[0] * du1[ii] +
                   z[1] * (-du1[1] * x_inv[ii - 2] +
                           v * du1[ii] * x_inv[ii - 2] * x_inv[ii - 2]));
      }
      Gdu[n + 2] = -(*dr);
      scs_float coeff0 = SCS(dot)(g0g1, du1, u1_dim) + g0g1[n + 2] * (*dr);
      scs_float coeff1 =
          SCS(dot)(g0g1 + u_dim, du1, u1_dim) + g0g1[2 * n + 5] * (*dr);
      memcpy(CCTdu, g0g1, u_dim * sizeof(*g0g1));
      SCS(scale_array)(CCTdu, coeff0, u_dim);
      SCS(add_scaled_array)(CCTdu, g0g1 + u_dim, u_dim, coeff1);
      CCTdu[n + 2] += (*dr);
      memcpy(residual, bnew, u_dim * sizeof(*bnew));

      // important to compute the residual as residual = bnew - (Gdu +
      // CCTdu) and not as residual = (bnew - Gdu) - CCTdu
      SCS(add_scaled_array)(Gdu, CCTdu, u_dim, 1.0);
      SCS(add_scaled_array)(residual, Gdu, u_dim, -1.0);
    }
  }

  // -------------------------------------------------------------------------
  //                 recover substituted variables
  // -------------------------------------------------------------------------
  memcpy(dz, rhs_reduced + u_dim, 3 * sizeof(*rhs_reduced));
  dz[0] -= w[0] * (SCS(dot)(g0g1, du1, u1_dim) + g0g1[n + 2] * (*dr));
  dz[1] -=
      w[1] * (SCS(dot)(g0g1 + u_dim, du1, u1_dim) + g0g1[2 * n + 5] * (*dr));
  dz[2] += du1[1];
  for (scs_int ii = 0; ii < 3; ++ii) {
    dz[ii] = -dz[ii] / (w[ii] * w[ii]);
    ds[ii] = w[ii] * (rhs2[ii] / lmbda[ii] - w[ii] * dz[ii]);
  }
}

/* A primal-dual interior point method for projecting the point (t0, v0, x0)
 * onto the logarithm cone. Workspace must be of size 22n + 122, where n is
 * the dimension of x0. The variable is u1 = [t, v, x], and epigraph
 * variable r
 */
scs_int log_cone_IPM(scs_float t0, scs_float v0, scs_float *x0, scs_float *u1,
                     scs_int n, scs_float *workspace, Newton_stats *stats,
                     scs_int variant) {
  scs_int m = 3;
  scs_int u_dim = n + 3;
  scs_int u1_dim = n + 2;

  // -----------------------------------------------------------------
  //                    scale (t0, v0, x0)
  // -----------------------------------------------------------------
  scs_float scale1 = SCS(norm_inf)(x0, n);
  scale1 = MAX(t0, scale1);
  scale1 = MAX(v0, scale1);
  scale1 = 1 / scale1;
  t0 = t0 * scale1;
  v0 = v0 * scale1;
  SCS(scale_array)(x0, scale1, n);

  // ----------------------------------------------------------------------
  //                       Parse workspace
  // ----------------------------------------------------------------------
  scs_float *x_inv = workspace;
  scs_float *z = workspace + n;
  scs_float *s = z + m;
  scs_float *f_u = s + m;
  scs_float *grad_f0 = f_u + m;
  scs_float *grad_f1 = grad_f0 + u_dim;
  scs_float *rx = grad_f1 + u_dim;
  scs_float *rznl = rx + u_dim;
  scs_float *w = rznl + m;
  scs_float *lmbda = w + m;
  scs_float *du1 = lmbda + m;
  scs_float *dz = du1 + u1_dim;
  scs_float *ds = dz + m;

  // for evaluating the new residuals in the line search
  scs_float *u1_new = ds + m;
  scs_float *z_new = u1_new + u1_dim;
  scs_float *s_new = z_new + m;
  scs_float *x_inv_new = s_new + m;
  scs_float *f_u_new = x_inv_new + n;
  scs_float *grad_f0_new = f_u_new + m;
  scs_float *grad_f1_new = grad_f0_new + u_dim;
  scs_float *rx_new = grad_f1_new + u_dim;
  scs_float *rznl_new = rx_new + u_dim;

  // for storing the state for the line search
  scs_float *u1_0 = rznl_new + m;
  scs_float *du1_0 = u1_0 + u1_dim;
  scs_float *z0 = du1_0 + u1_dim;
  scs_float *dz0 = z0 + m;
  scs_float *s0 = dz0 + m;
  scs_float *ds0 = s0 + m;
  scs_float r0 = 0.0;
  scs_float dr0 = 0.0;

  KKT_solver_workspace KKT_work;
  KKT_work.g0g1 = grad_f0;
  KKT_work.temp1 = ds0 + m;
  KKT_work.temp2 = KKT_work.temp1 + u_dim;
  KKT_work.GinvC = KKT_work.temp2 + u_dim;
  KKT_work.R = KKT_work.GinvC + 3 * u_dim;
  KKT_work.GinvRes = KKT_work.R + 9;
  KKT_work.CTGinvRes = KKT_work.GinvRes + u_dim;
  KKT_work.RinvCTGinvRes = KKT_work.CTGinvRes + 3;
  KKT_work.residual = KKT_work.RinvCTGinvRes + 3;
  KKT_work.rhs_reduced = KKT_work.residual + u_dim;
  KKT_work.bnew = KKT_work.rhs_reduced + u_dim + 3;
  KKT_work.CCTdu = KKT_work.bnew + u_dim;

  // ---------------------------------------------------------------------
  // initialize primal variable u1 = [t, v, x], epigraph variable r and
  // primal slack variable s and Lagrange multipliers z
  // ---------------------------------------------------------------------
  for (scs_int i = 0; i < n + 2; ++i) {
    u1[i] = 1.0;
  }

  scs_float r = 0.0;
  scs_float dr = 0.0;
  scs_float rnew = 0.0;

  for (scs_int i = 0; i < m; ++i) {
    z[i] = 1.0;
    s[i] = 1.0;
  }

  scs_int relaxed_iters = 0;
  size_t iter = 0;
  scs_float theta1, theta2, theta3;
  scs_float pres0, dres0;
  scs_float phi0 = 0.0, dphi0 = 0.0, step_size0 = 0.0;
  // scs_int small_consecutive_steps_counter = 0;

#ifdef SPECTRAL_DEBUG
  printf("%-3s%-15s%-15s%-15s%-10s%-10s\n", "", "gap", "pres", "dres", "sig",
         "step");
#endif

  for (iter = 0; iter < MAX_ITER_IPM; ++iter) {
    scs_float gap = z[0] * s[0] + z[1] * s[1] + z[2] * s[2];
    scs_float mu = gap / m;

    // --------------------------------------------------------
    //                   evaluate oracle
    // --------------------------------------------------------
    for (scs_int i = 0; i < n; i++) {
      x_inv[i] = 1 / u1[i + 2];
    }
    f_oracle(u1, r, t0, v0, x0, x_inv, n, f_u, grad_f0, grad_f1);

    // --------------------------------------------------------
    //    compute residuals and check termination criteria
    // --------------------------------------------------------
    memcpy(rx, grad_f0, u_dim * sizeof(*grad_f0));
    SCS(scale_array)(rx, z[0], u_dim);
    SCS(add_scaled_array)(rx, grad_f1, u_dim, z[1]);
    rx[1] -= z[2];
    rx[n + 2] += 1;
    rznl[0] = f_u[0] + s[0];
    rznl[1] = f_u[1] + s[1];
    rznl[2] = f_u[2] + s[2];

    scs_float dres = SCS(norm_2)(rx, u_dim);
    scs_float pres =
        sqrt(rznl[0] * rznl[0] + rznl[1] * rznl[1] + rznl[2] * rznl[2]);
    scs_float norm_rx = dres;
    scs_float norm_rznl = pres;

    if (iter == 0) {
      pres0 = MAX(pres, 1.0);
      dres0 = MAX(dres, 1.0);
      theta1 = 1.0 / gap;
      theta2 = 1.0 / dres0;
      theta3 = 1.0 / pres0;
    }

    pres = pres / pres0;
    dres = dres / dres0;
    scs_float relgap = gap / MAX(r, 1.0);

    if (dres < FEASTOL_IPM && pres < FEASTOL_IPM &&
        (gap < ABSTOL_IPM || relgap <= RELTOL_IPM)) {
#ifdef SPECTRAL_DEBUG
      printf("optimal solution found: \n");
      printf("gap / pres / dres: %.7e, %.7e, %.7e \n", gap, pres, dres);
#endif
      break;
    }

    // ----------------------------------------------------
    //      compute scaling matrix and scaling point
    // ----------------------------------------------------
    for (scs_int i = 0; i < m; i++) {
      w[i] = sqrt(s[i] / z[i]);
      lmbda[i] = sqrt(s[i] * z[i]);
    }

    // -----------------------------------------------------------------
    //      precomputations for KKT system
    // -----------------------------------------------------------------
    scs_float scale2 = 1 / w[0];
    SCS(scale_array)(grad_f0, scale2, u_dim);
    scale2 = 1 / w[1];
    SCS(scale_array)(grad_f1, scale2, u_dim);
    KKT_precompute(u1[1], x_inv, z, w, n, &KKT_work);

    scs_float rhs2_aff[] = {-lmbda[0] * lmbda[0], -lmbda[1] * lmbda[1],
                            -lmbda[2] * lmbda[2]};

    // upper bound on loop is chosen so it also iterates over rznl
    for (scs_int i = 0; i < u_dim + 3; ++i) {
      rx[i] = -rx[i];
    }
    scs_float *rhs1_aff = rx;

    scs_float sigma = 0.0;
    scs_float phi = 0.0;
    scs_float dphi = 0.0;
    scs_float step_size = 0.0;
    phi = theta1 * gap + theta2 * norm_rx + theta3 * norm_rznl;
    dphi = -theta1 * (1 - sigma) * gap - theta2 * norm_rx - theta3 * norm_rznl;

    for (scs_int i = 0; i < 2; ++i) {
      // ------------------------------------------------------------
      //  For i = 0 we compute the affine-scaling direction.
      //  For i = 1 we compute the actual search direction.
      // ------------------------------------------------------------
      if (i == 1 && variant == 0) {
        SCS(scale_array)(rhs1_aff, 1 - sigma, u_dim + 3);
        rhs2_aff[0] += (sigma * mu - ds[0] * dz[0]);
        rhs2_aff[1] += (sigma * mu - ds[1] * dz[1]);
        rhs2_aff[2] += (sigma * mu - ds[2] * dz[2]);
      }

      KKT_solve(z, w, x_inv, u1[1], lmbda, n, rhs1_aff, rhs2_aff, &KKT_work,
                du1, &dr, dz, ds);

      // --------------------------------------------------------------
      //  For i = 0 we determine the centering parameter.
      //  For i = 1 we do a nonmonotone line search.
      // --------------------------------------------------------------
      step_size = find_max_step_size(u1, z, s, du1, dz, ds, n);

      bool backtrack = true;

      while (backtrack) {
        // u_new = u + step * du, z_new = z + step * dz, s_new = s + step
        memcpy(u1_new, u1, u1_dim * sizeof(*u1));
        SCS(add_scaled_array)(u1_new, du1, u1_dim, step_size);
        rnew = r + step_size * dr;
        memcpy(z_new, z, 3 * sizeof(*z));
        SCS(add_scaled_array)(z_new, dz, m, step_size);
        memcpy(s_new, s, 3 * sizeof(*s));
        SCS(add_scaled_array)(s_new, ds, m, step_size);

        // evaluate oracle in u_new
        for (scs_int i = 0; i < n; i++) {
          x_inv_new[i] = 1 / u1_new[i + 2];
        }
        f_oracle(u1_new, rnew, t0, v0, x0, x_inv_new, n, f_u_new, grad_f0_new,
                 grad_f1_new);

        // compute residuals and merit function
        memcpy(rx_new, grad_f0_new, u_dim * sizeof(*grad_f0_new));
        SCS(scale_array)(rx_new, z_new[0], u_dim);
        SCS(add_scaled_array)(rx_new, grad_f1_new, u_dim, z_new[1]);
        rx_new[1] -= z_new[2];
        rx_new[n + 2] += 1;
        rznl_new[0] = f_u_new[0] + s_new[0];
        rznl_new[1] = f_u_new[1] + s_new[1];
        rznl_new[2] = f_u_new[2] + s_new[2];
        scs_float gap_new = SCS(dot)(z_new, s_new, m);
        scs_float phi_new = theta1 * gap_new +
                            theta2 * SCS(norm_2)(rx_new, u_dim) +
                            theta3 * SCS(norm_2)(rznl_new, m);
        // ----------------------------------------------------------
        //   For i == 0 we determine the centering parameter
        // ----------------------------------------------------------
        if (i == 0) {
          if (phi_new <= (1 - ALPHA_IPM * step_size) * phi) {
            backtrack = false;
            sigma = (gap_new / gap);
            if (sigma < 1) {
              sigma *= (sigma * sigma);
            }
          } else {
            step_size *= BETA_IPM;
          }
        }
        // ------------------------------------------------------------
        //  For i == 1 we do a nonmonotone line search with the actual
        //  search direction
        // ------------------------------------------------------------
        else {
          if (relaxed_iters == -1 || MAX_RELAXED_ITERS == 0) {
            if (phi_new <= phi + ALPHA_IPM * step_size * dphi) {
              // relaxed_iters = 0;
              backtrack = false;
            } else {
              step_size *= BETA_IPM;
            }
          } else if (relaxed_iters == 0) {
            if (phi_new <= phi + ALPHA_IPM * step_size * dphi) {
              relaxed_iters = 0;
            } else {
              relaxed_iters = 1;
              memcpy(u1_0, u1, u1_dim * sizeof(*u1));
              r0 = r;
              memcpy(z0, z, 3 * sizeof(*z));
              memcpy(s0, s, 3 * sizeof(*s));
              memcpy(du1_0, du1, u1_dim * sizeof(*du1));
              dr0 = dr;
              memcpy(dz0, dz, 3 * sizeof(*dz));
              memcpy(ds0, ds, 3 * sizeof(*ds));
              phi0 = phi;
              dphi0 = dphi;
              step_size0 = step_size;
            }

            backtrack = false;
          } else if (relaxed_iters == MAX_RELAXED_ITERS) {
            if (phi_new <= phi0 + ALPHA_IPM * step_size0 * dphi0) {
              // relaxed_iters = 0;
              backtrack = false;
            } else {
              relaxed_iters = -1;
              memcpy(u1, u1_0, u1_dim * sizeof(*u1));
              r = r0;
              memcpy(z, z0, 3 * sizeof(*z));
              memcpy(s, s0, 3 * sizeof(*s));
              memcpy(du1, du1_0, u1_dim * sizeof(*du1));
              dr = dr0;
              memcpy(dz, dz0, 3 * sizeof(*dz));
              memcpy(ds, ds0, 3 * sizeof(*ds));
              phi = phi0;
              dphi = dphi0;
              step_size = step_size0;
            }
          } else if (relaxed_iters < MAX_RELAXED_ITERS) {
            if (phi_new <= phi0 + ALPHA_IPM * step_size0 * dphi0) {
              relaxed_iters = 0;
            } else {
              relaxed_iters += 1;
            }

            backtrack = false;
          }
        }
      }
    }

    // update iterates
    memcpy(u1, u1_new, u1_dim * sizeof(*u1));
    r = rnew;
    memcpy(z, z_new, 3 * sizeof(*z));
    memcpy(s, s_new, 3 * sizeof(*s));
#ifdef SPECTRAL_DEBUG
    printf("%ld: %.7e, %.7e, %.7e, %f, %.3f \n", iter, gap, pres, dres, sigma,
           step_size);
#endif
  }

  // unscale solution
  scale1 = 1 / scale1;
  SCS(scale_array)(x0, scale1, n);
  SCS(scale_array)(u1, scale1, u1_dim);
  stats->iter = iter;
  return 0;
}
