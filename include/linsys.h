#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

/* This is the API that any new linear system solver must implement */

/* Struct containing linear system workspace. Implemented by linear solver. */
/* This typedef is in scs.h */
/* typedef struct SCS_LIN_SYS_WORK ScsLinSysWork; */

/**
 * Initialize `ScsLinSysWork` structure and perform any necessary preprocessing.
 *
 *  @param  A          `A` data matrix, `m x n`.
 *  @param  P          `P` data matrix, `n x n`.
 *  @param  diag_r     `R > 0` diagonal entries of length `m + n`.
 *  @return            Linear system solver workspace.
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
 * Solves the linear system as required by SCS at each iteration:
 * \f[
 *    \begin{bmatrix}
 *    (R_x + P) & A^\top \\
 *     A   &  -R_y \\
 *    \end{bmatrix} x = b
 *  \f]
 *
 *  for `x`, where `diag(R_x, R_y) = R`. Overwrites `b` with result.
 *
 *  @param  w    Linear system private workspace.
 *  @param  b    Right hand side, contains solution at the end.
 *  @param  s    Contains warm-start (may be NULL).
 *  @param  tol  Tolerance required for the system solve.
 *  @return status != 0 indicates failure.
 *
 */
scs_int scs_solve_lin_sys(ScsLinSysWork *w, scs_float *b, const scs_float *s,
                          scs_float tol);
/**
 *  Update the linsys workspace when `R` is changed. For example, a
 *  direct method for solving the linear system might need to update the
 *  factorization of the matrix.
 *
 *  @param  w             Linear system private workspace.
 *  @param  new_diag_r    Updated `diag_r`, diagonal entries of R.
 *
 */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *w, const scs_float *new_diag_r);

/*
 * Optional: products with the linear solver's own copy of the matrix.
 *
 * Direct solvers refine each iteration's solve in correction form while the
 * solve's error limits SCS (see "Accuracy of each iteration's solve" in
 * glbopts.h): from the previous iteration's solution x0 they solve
 * K d = b - K x0 and take x = x0 + d, where
 * \f$K = \begin{bmatrix} R_x + P & A^\top \\ A & -R_y \end{bmatrix}\f$
 * with the current `diag_r`. That costs one multiply by K per iteration, and
 * SCS checks the residual b - K x of a solve every 25 iterations. By default
 * SCS computes these products itself, from the problem data. A solver that
 * holds K anyway can compute them with its own copy instead, for instance
 * with a threaded or GPU sparse matrix-vector product, by implementing the
 * three functions below. SCS makes every decision (when to check, when to
 * use the correction form, the thresholds); the solver only computes. All
 * vectors are host vectors of length n + m.
 *
 * Contract:
 *  - Whether a workspace provides the products is fixed when
 *    `scs_init_lin_sys_work` returns, by a rule stated up front (for instance
 *    the size of K), never by whether a resource happened to be available.
 *    A workspace that does not provide them returns
 *    `SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED` from every `scs_lin_sys_residual`
 *    call and does nothing in the other two functions.
 *  - A workspace that should provide them but cannot set them up (out of
 *    memory, a library error, K too large for the product's index type) makes
 *    `scs_init_lin_sys_work` fail, with a message saying why; it does not
 *    fall back to "not provided".
 *  - Any error in a product returns `SCS_LIN_SYS_PRODUCTS_FAILED`.
 *  - SCS calls `scs_lin_sys_residual` once at setup, right after
 *    `scs_init_lin_sys_work`, with zero vectors: "not provided" makes SCS
 *    compute every product itself, a failure (any negative value other than
 *    "not provided") makes SCS's setup fail, and a norm (>= 0) means the
 *    solver provides them. During a solve every negative return is a failure
 *    of the linear solver, and the solve stops.
 */

/** Returned by `scs_lin_sys_residual`: this workspace does not compute the
 * products, SCS computes them itself. */
#define SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED (-1.)
/** Returned by `scs_lin_sys_residual`: the product failed. */
#define SCS_LIN_SYS_PRODUCTS_FAILED (-2.)

/**
 * Computes \f$r = b - K x\f$ with the solver's own copy of K (as updated by
 * the last `scs_update_lin_sys_diag_r`) and returns \f$\|r\|_\infty\f$. The
 * solver keeps r for a following `scs_lin_sys_correct`. If `x` is `SCS_NULL`,
 * x is the vector stored by the last `scs_lin_sys_keep`. It must not change
 * the solver's copy of the solution of its latest solve (see
 * `scs_lin_sys_keep`).
 *
 *  @param  w    Linear system private workspace.
 *  @param  b    Right-hand side.
 *  @param  x    Point at which the residual is evaluated, or `SCS_NULL`.
 *  @return      \f$\|b - K x\|_\infty\f$,
 *               `SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED` or
 *               `SCS_LIN_SYS_PRODUCTS_FAILED` (see the contract above).
 */
scs_float scs_lin_sys_residual(ScsLinSysWork *w, const scs_float *b,
                               const scs_float *x);

/**
 * Stores `x` as the guess x0 of the next correction-form solve. `x` is always
 * the solution returned by the latest `scs_solve_lin_sys` or
 * `scs_lin_sys_correct` call, unchanged, so a solver that keeps that solution
 * on its side (a GPU solver) may store its own copy and not read `x`.
 * Between that solve and this call SCS may call `scs_lin_sys_residual` (with
 * that solution or another vector as `x`, or with `SCS_NULL`), but no other
 * solve: the residual must therefore leave the solver's copy of the solution
 * as it is. Called only if the solver provides the products.
 *
 *  @param  w    Linear system private workspace.
 *  @param  x    The solution of the latest solve.
 */
void scs_lin_sys_keep(ScsLinSysWork *w, const scs_float *x);

/**
 * Solves in correction form: \f$x = x_0 + K^{-1} r\f$, where r is the
 * residual of the latest `scs_lin_sys_residual` call, which SCS always makes
 * with `x` equal to `SCS_NULL` (so \f$r = b - K x_0\f$) right before. The
 * result counts as the latest solve for `scs_lin_sys_keep`. Called only if
 * the solver provides the products.
 *
 *  @param  w    Linear system private workspace.
 *  @param  x    Output, the solution.
 *  @return      status != 0 indicates failure.
 */
scs_int scs_lin_sys_correct(ScsLinSysWork *w, scs_float *x);

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
