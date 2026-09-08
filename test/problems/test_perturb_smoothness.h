#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Perturbation smoothness: is the map b -> xstar(b) numerically
 * differentiable?
 *
 * Downstream users differentiate through SCS (diffcp, cvxpylayers). Their
 * gradcheck compares an analytic derivative against a finite difference, so
 * they need (xstar(b + h e_j) - xstar(b)) / h to be a usable estimate of the
 * partial derivative of xstar with respect to b_j. That holds only if two
 * nearby solves make *correlated* errors: each is accurate to ~eps, and h is
 * small, so uncorrelated errors of size eps enter the quotient magnified by
 * 1/h.
 *
 * Several things could decorrelate them. Since 3.3.0 the Ruiz equilibration
 * is a function of b and c as well as A and P (see the tau-column seeding of
 * D in compute_l2_mats), so a perturbed solve is scaled slightly differently;
 * adaptive_diag_scale likewise refines the metric from observed residual
 * profiles. Both are solve-path dependent in a way the pre-3.3.0 code was not.
 *
 * Measured effect on the two problems below: none that matters. This test is
 * therefore a *regression guard*, not a reproduction of a known defect -- it
 * pins the property that downstream differentiation depends on, so that a
 * future change which does break it fails here rather than in someone else's
 * gradcheck.
 *
 * Two problems, two ways of measuring, because neither alone covers the code
 * that 3.3.0 changed:
 *
 * 1. Equality-constrained QP (zero cone only), exact Jacobian.
 *
 *      minimize 0.5 ||x||^2  s.t.  x0 + x1 = b0,  x1 + x2 = b1
 *
 *    With Minv = inv(A A^T) = (1/3) [[2, -1], [-1, 2]], the solution is
 *    xstar(b) = A^T Minv b, hence d xstar / d b0 = (1/3) [2, 1, -1] exactly,
 *    for every b. The map is affine, so the forward difference carries no
 *    truncation error and any deviation is pure numerical noise.
 *    b = [1, 1000] is badly scaled on purpose, so the equilibration has real
 *    work to do and its b-dependence is not masked.
 *
 * 2. SDP with a PSD block, self-consistency -- compiled only under
 *    USE_LAPACK, since the PSD projection needs BLAS/LAPACK and SCS refuses
 *    the cone outright without it. The zero cone never engages the
 *    per-row metric refinement or the PSD block metric, which is where the
 *    3.3.0 changes live, so the QP alone would not exercise them.
 *
 *      minimize <C, X>  s.t.  trace(X) = b0,  X PSD  (2x2)
 *
 *    The solution is X = b0 * v v^T with v the min-eigenvector of C, so xstar
 *    is again exactly linear in b0. Here we compare the forward difference at
 *    h against the one at h/2: for an affine map they agree exactly, so their
 *    discrepancy needs no analytic Jacobian and no assumption about the svec
 *    convention.
 *
 * Thresholds are absolute and loose. Measured across the direct, indirect,
 * dense and accelerate backends, the worst normalize=1 figures are 2.3e-06
 * (QP) and 8.8e-08 (SDP), both on the indirect solver at h=1e-06; the bound
 * of 1e-04 therefore keeps roughly a factor of 40 in hand. These errors are
 * dominated by the eps/h floor inherent to finite differencing rather than by
 * anything SCS chooses, and a genuine loss of smoothness shows up orders of
 * magnitude above them. Asserting a ratio against normalize=0 instead was
 * tried and rejected: the denominator sits at that same floor, so the ratio
 * swings between 0.3 and 11.5 across regimes where the absolute error barely
 * moves.
 */

#define PS_N 3

/* --- 1. equality-constrained QP ------------------------------------------ */

/* exact d xstar / d b0, independent of b */
static const scs_float ps_jac_b0[PS_N] = {2.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0};

/* Solve with b0 shifted by db0. Fresh data each call: SCS normalizes A, P, b
 * and c in place. */
static scs_int ps_solve_qp(scs_float db0, scs_int normalize,
                           scs_float *x_out) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag, i;

  scs_float Ax[] = {1.0, 1.0, 1.0, 1.0};
  scs_int Ai[] = {0, 0, 1, 1};
  scs_int Ap[] = {0, 1, 3, 4};
  scs_float Px[] = {1.0, 1.0, 1.0};
  scs_int Pi[] = {0, 1, 2};
  scs_int Pp[] = {0, 1, 2, 3};
  scs_float b[] = {1.0, 1000.0};
  scs_float c[] = {0.0, 0.0, 0.0};

  b[0] += db0;

  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  d->m = 2;
  d->n = PS_N;
  d->b = b;
  d->c = c;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = 2;
  d->A->n = PS_N;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;
  d->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->P->m = PS_N;
  d->P->n = PS_N;
  d->P->x = Px;
  d->P->i = Pi;
  d->P->p = Pp;
  k->z = 2;

  scs_set_default_settings(stgs);
  stgs->normalize = normalize;
/* the indirect (CG) backend does not reach 1e-12 here; errors
   * saturate by 1e-8 anyway, so 1e-9 costs no resolution */
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;
  stgs->max_iters = 200000;
  stgs->verbose = 0;

  exitflag = scs(d, k, stgs, sol, &info);
  for (i = 0; i < PS_N; ++i) {
    x_out[i] = sol->x[i];
  }

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d->P);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return exitflag;
}

/* max-norm error of the forward difference against the exact Jacobian */
static scs_float ps_qp_fd_err(scs_float h, scs_int normalize) {
  scs_float x0[PS_N], x1[PS_N], err = 0.0, fd;
  scs_int i;

  if (ps_solve_qp(0.0, normalize, x0) != SCS_SOLVED) {
    return -1.0;
  }
  if (ps_solve_qp(h, normalize, x1) != SCS_SOLVED) {
    return -1.0;
  }
  for (i = 0; i < PS_N; ++i) {
    fd = (x1[i] - x0[i]) / h;
    if (ABS(fd - ps_jac_b0[i]) > err) {
      err = ABS(fd - ps_jac_b0[i]);
    }
  }
  return err;
}

/* --- 2. SDP with a PSD block --------------------------------------------- */

/* PSD projection needs BLAS/LAPACK; without it SCS refuses the cone outright
 * ("FATAL: SDP/Complex SDP requires BLAS/LAPACK"), so this half is compiled
 * only when it can run. The QP half above needs no LAPACK and always runs. */
#if defined(USE_LAPACK)

/* Rows: 0 = trace equality (z=1), 1..3 = PSD block (s=[2]). */
static scs_int ps_solve_sdp(scs_float db0, scs_int normalize,
                            scs_float *x_out) {
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag, i;

  scs_float Ax[] = {1.0, -1.0, -1.0, 1.0, -1.0};
  scs_int Ai[] = {0, 1, 2, 0, 3};
  scs_int Ap[] = {0, 2, 3, 5};
  scs_float b[] = {1.0, 0.0, 0.0, 0.0};
  /* svec(C) for C = [[1, 0.5], [0.5, 2]] */
  scs_float c[] = {1.0, 0.5 * 1.41421356237309505, 2.0};
  scs_int sd[] = {2};

  b[0] += db0;

  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  d->m = 4;
  d->n = PS_N;
  d->b = b;
  d->c = c;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = 4;
  d->A->n = PS_N;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;
  k->z = 1;
  k->s = sd;
  k->ssize = 1;

  scs_set_default_settings(stgs);
  stgs->normalize = normalize;
  /* 1e-12 is not attainable on this SDP; 1e-9 is, comfortably */
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;
  stgs->max_iters = 200000;
  stgs->verbose = 0;

  exitflag = scs(d, k, stgs, sol, &info);
  for (i = 0; i < PS_N; ++i) {
    x_out[i] = sol->x[i];
  }

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return exitflag;
}

/* |fd(h) - fd(h/2)|_inf; exactly zero for an affine solution map */
static scs_float ps_sdp_fd_incons(scs_float h, scs_int normalize) {
  scs_float x0[PS_N], xa[PS_N], xb[PS_N], err = 0.0, qa, qb;
  scs_int i;

  if (ps_solve_sdp(0.0, normalize, x0) != SCS_SOLVED) {
    return -1.0;
  }
  if (ps_solve_sdp(h, normalize, xa) != SCS_SOLVED) {
    return -1.0;
  }
  if (ps_solve_sdp(h / 2, normalize, xb) != SCS_SOLVED) {
    return -1.0;
  }
  for (i = 0; i < PS_N; ++i) {
    qa = (xa[i] - x0[i]) / h;
    qb = (xb[i] - x0[i]) / (h / 2);
    if (ABS(qa - qb) > err) {
      err = ABS(qa - qb);
    }
  }
  return err;
}

#endif /* USE_LAPACK */

/* --- the test ------------------------------------------------------------ */

static const char *test_perturb_smoothness(void) {
  const scs_float hs[] = {1e-4, 1e-6};
  const scs_int nh = 2;
  scs_int j;

  for (j = 0; j < nh; ++j) {
    scs_float h = hs[j];
    scs_float qp_off = ps_qp_fd_err(h, 0);
    scs_float qp_on = ps_qp_fd_err(h, 1);

    mu_assert("perturb_smoothness: QP solve failed", qp_off >= 0 && qp_on >= 0);

    scs_printf("perturb_smoothness: h=%.0e  QP fd_err norm0/norm1 %.2e/%.2e\n",
               h, qp_off, qp_on);

    /* normalize=1 is the shipped default and the one downstream code hits */
    mu_assert("perturb_smoothness: QP finite difference lost accuracy -- the "
              "solution map is no longer smooth enough to differentiate",
              qp_on < 1e-4);

#if defined(USE_LAPACK)
    {
      scs_float sdp_off = ps_sdp_fd_incons(h, 0);
      scs_float sdp_on = ps_sdp_fd_incons(h, 1);

      mu_assert("perturb_smoothness: SDP solve failed",
                sdp_off >= 0 && sdp_on >= 0);

      scs_printf("perturb_smoothness: h=%.0e  SDP fd_incons norm0/norm1 "
                 "%.2e/%.2e\n",
                 h, sdp_off, sdp_on);

      mu_assert("perturb_smoothness: SDP finite difference is inconsistent "
                "between step sizes -- PSD path lost smoothness",
                sdp_on < 1e-4);
    }
#endif
  }
  return 0;
}
