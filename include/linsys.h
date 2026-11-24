#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

/* This is the API that any new linear system solver must implement. See
 * "Implementing a new linear solver" in the documentation: the solvers in
 * linsys/ implement scs_solve_lin_sys as one call to a shared routine
 * (SCS(lin_sys_shared_solve), linsys/scs_matrix.h) around their plain solve,
 * and a new solver can do the same. */

/* Struct containing linear system workspace. Implemented by linear solver. */
/* This typedef is in scs.h */
/* typedef struct SCS_LIN_SYS_WORK ScsLinSysWork; */

/**
 * Initialize `ScsLinSysWork` structure and perform any necessary preprocessing
 * (for instance factorize K, see `scs_solve_lin_sys`).
 *
 * The solver may keep the three pointers and read the arrays in any later
 * call: they stay valid, and A and P unchanged, until
 * `scs_free_lin_sys_work`, and SCS writes every new R into this same
 * `diag_r` array before it calls `scs_update_lin_sys_diag_r`.
 *
 *  @param  A          `A` data matrix, `m x n`.
 *  @param  P          `P` data matrix, `n x n`, upper triangle, or `SCS_NULL`.
 *  @param  diag_r     `R > 0` diagonal entries of length `m + n`.
 *  @return            Linear system solver workspace, or `SCS_NULL` on
 *                     failure (with a message saying why).
 *
 */
ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r);

/**
 * Frees `ScsLinSysWork` structure and associated allocated memory.
 *
 *  @param  w    Linear system private workspace.
 */
void scs_free_lin_sys_work(ScsLinSysWork *w);

/**
 * Solves the linear system as required by SCS:
 * \f[
 *    K x = b, \qquad
 *    K = \begin{bmatrix}
 *    (R_x + P) & A^\top \\
 *     A   &  -R_y \\
 *    \end{bmatrix}
 *  \f]
 *
 *  for `x`, where `diag(R_x, R_y) = R` is the one given to the last
 *  `scs_update_lin_sys_diag_r` (or to `scs_init_lin_sys_work`), and
 *  overwrites `b` with `x`. With `max_refine_steps` 0 and `res` `SCS_NULL`
 *  this is the solve of SCS 3.3 for a CG solver (conjugate gradients from
 *  `s` to `tol`), and for a direct solver given no `s` (a plain solve).
 *
 *  Direct solvers (a factorization: QDLDL, dense, Apple Accelerate, MKL
 *  Pardiso, cuDSS):
 *  1. Warm start. Given `s`, the solver computes \f$r = b - K s\f$ and, if
 *     \f$\|r\|_\infty < \|b\|_\infty\f$, returns the correction form
 *     \f$x = s + K^{-1} r\f$, whose error is about that of a plain solve
 *     times \f$\|r\|_\infty / \|b\|_\infty\f$, for one multiply by K more.
 *     Otherwise (s is no better a guess than zero, for instance after R
 *     changed) it solves plainly, \f$x = K^{-1} b\f$.
 *  2. Refinement, if `max_refine_steps` > 0: while
 *     \f$\|b - K x\|_\infty\f$ > `tol` and fewer than `max_refine_steps`
 *     steps were taken, the candidate \f$c = x + K^{-1}(b - K x)\f$
 *     replaces x if its residual is strictly smaller; if not, x stays and
 *     refinement stops. It also stops after a kept step that did not at
 *     least halve the residual. Each step costs a solve and a multiply by K,
 *     the first one a multiply more (the residual of x).
 *
 *  CG solvers (conjugate gradients on the reduced system, CPU and GPU):
 *  1. CG from `s` (its first n entries; y is computed from x), or from
 *     zero, until the residual CG tracks is below `tol`.
 *  2. Refinement, if `max_refine_steps` > 0: as for the direct solvers, with
 *     \f$K^{-1}(b - K x)\f$ a CG solve from zero to `tol`. The residual is
 *     computed in the full space with the multiplies below: the one CG
 *     tracks, on the reduced system, can drift from it (in single precision
 *     or from a warm start far from x) and has a higher rounding floor, and
 *     y is computed from x, so restarting CG on b from x would remove only
 *     part of either. Each step costs CG iterations and a multiply by K.
 *
 *  The multiplies by K (warm start, refinement, `res`) use the solver's own
 *  copy of K where it has one, by a rule fixed when `scs_init_lin_sys_work`
 *  returns and stated in the solver's documentation (MKL Pardiso; cuDSS from
 *  `CUDSS_MIN_PRODUCT_NNZ` nonzeros), else SCS's products with the A, P and
 *  diag_r given to `scs_init_lin_sys_work`. Every decision compares
 *  residuals, never timings: a call is deterministic whenever the solver's
 *  plain solve and products are.
 *
 *  @param  w    Linear system private workspace.
 *  @param  b    Right-hand side, n + m entries; overwritten with the
 *               solution.
 *  @param  s    Warm start, a guess of the solution (n + m entries; CG reads
 *               the first n), or `SCS_NULL`. Not changed; must not alias `b`.
 *  @param  tol  The accuracy asked for, as \f$\|b - K x\|_\infty\f$. CG
 *               solvers stop at it, also below 1e-12 (smaller is more
 *               accurate and costs more CG iterations), and their refinement
 *               steps solve to it; they need `tol` > 0, and return x = 0
 *               without iterating when
 *               \f$\|b\|_\infty \le \min(\mathrm{tol}, 10^{-12})\f$.
 *               Direct solvers refine toward it, and ignore it when
 *               `max_refine_steps` is 0. Below the rounding level of the
 *               residual (at most a few machine epsilons times
 *               \f$\| |b| + |K|\,|x| \|_\infty\f$, absolute values entry by
 *               entry, and typically a third of one, since the rounding
 *               errors of a row's terms partly cancel) a step only measures
 *               rounding, and the halving test stops it.
 *  @param  max_refine_steps  At most this many refinement steps; 0 for none.
 *  @param  res  If not `SCS_NULL`, set to \f$\|b - K x\|_\infty\f$ for the
 *               returned x and the given b, with the multiplies above.
 *               Cost: a copy of b and a multiply by K, unless refinement
 *               computed it. Asking for it never changes x. A NaN entry of
 *               the residual may be missed by the max (as by SCS's
 *               `norm_inf`), so a finite `*res` does not show that x is
 *               finite.
 *  @return 0 on success; nonzero if a solve or a multiply failed (b and
 *          `*res` are then unspecified).
 *
 */
scs_int scs_solve_lin_sys(ScsLinSysWork *w, scs_float *b, const scs_float *s,
                          scs_float tol, scs_int max_refine_steps,
                          scs_float *res);
/**
 *  Update the linsys workspace when `R` is changed. For example, a
 *  direct method for solving the linear system might need to update the
 *  factorization of the matrix. SCS passes the `diag_r` array given to
 *  `scs_init_lin_sys_work`, with the new R written into it.
 *
 *  @param  w             Linear system private workspace.
 *  @param  new_diag_r    Updated `diag_r`, diagonal entries of R.
 *  @return 0 on success; any other value, negative or positive, is a failure
 *          (for instance a factorization that failed), after the solver has
 *          printed why. SCS then goes back to the previous R: it writes that
 *          R into the same array and calls this function again, so a failed
 *          update must leave the workspace able to take another one. If that
 *          call fails too, SCS stops with `SCS_FAILED`.
 *
 */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *w, const scs_float *new_diag_r);

/**
 * Name of the linear solver.
 *
 * @return name of method.
 */
const char *scs_get_lin_sys_method(void);

#ifdef __cplusplus
}
#endif

#endif
