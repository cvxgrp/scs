/*
 * Sparse matrix operations: normalization/equilibration, copying, validation,
 * and matrix-vector products (A*x, A'*x, P*x) for CSC-format matrices.
 */

#ifndef SCS_MATRIX_H_GUARD
#define SCS_MATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "linsys.h"
#include "scs_work.h"

/* Normalization routines, used if d->NORMALIZE is true */
/* normalizes A matrix, sets scal->E and scal->D diagonal scaling matrices,
 * A -> D*A*E. D and E must be all positive entries, D must satisfy cone
 * boundaries */
ScsScaling *SCS(normalize_a_p)(ScsMatrix *P, ScsMatrix *A, const scs_float *b,
                               const scs_float *c, ScsConeWork *cone);

/* to free the memory allocated in a ScsMatrix (called on A and P at finish) */
void SCS(free_scs_matrix)(ScsMatrix *A);

/* copies A (instead of in-place normalization), returns 0 for failure,
 * allocates memory for dstp	*/
scs_int SCS(copy_matrix)(ScsMatrix **dstp, const ScsMatrix *src);

/**
 * Validate the linear system inputs, returns < 0 if not valid inputs.
 *
 *  @param  A    A data matrix
 *  @param  P    P data matrix
 *  @return status < 0 indicates failure
 */
scs_int SCS(validate_lin_sys)(const ScsMatrix *A, const ScsMatrix *P);

/**
 * Forms y += A^T * x
 *
 *  @param  A    A data matrix
 *  @param  x    Input
 *  @param  y    Output
 */
void SCS(accum_by_atrans)(const ScsMatrix *A, const scs_float *x, scs_float *y);

/**
 * Forms y += A * x
 *
 *  @param  A           Data matrix
 *  @param  x           Input
 *  @param  y           Output
 */
void SCS(accum_by_a)(const ScsMatrix *A, const scs_float *x, scs_float *y);

/**
 * Forms y += P * x
 *
 *  @param  P    P data matrix
 *  @param  x    Input
 *  @param  y    Output
 */
void SCS(accum_by_p)(const ScsMatrix *P, const scs_float *x, scs_float *y);

/* ------------------------------------------------------------------------
 * The shared implementation of scs_solve_lin_sys (linsys.h).
 *
 * Every linear solver in linsys/ implements scs_solve_lin_sys as one call to
 * SCS(lin_sys_shared_solve), with an ScsLinSysShared that
 * SCS(lin_sys_shared_init) sets up around the solver's plain solve. The
 * routine adds what linsys.h asks of the solve beyond the plain solve: the
 * warm start (direct solvers: the correction form and its comparison with
 * zero; CG: passed on to CG), the refinement steps (direct) or restarts
 * (CG), and the residual. Its multiplies by K are SCS's own products with
 * the A, P and diag_r given to scs_init_lin_sys_work, unless the solver
 * passes a residual op with its own copy of K (MKL Pardiso), and with it,
 * optionally, a correct op that keeps the residual where the solver keeps
 * its vectors (cuDSS, on the GPU). This header is not installed: the routine
 * is for solvers built with SCS's sources (linsys/scs_matrix.c).
 * ---------------------------------------------------------------------- */

/* LIN_SYS_REFINE_CONTRACTION: refinement steps (direct solvers) and restarts
 * (CG) inside scs_solve_lin_sys stop after a kept one that did not reduce
 * ||b - K x||_inf to at most this fraction of the previous residual. That is
 * LAPACK's stagnation test (xGERFS): a step that does not halve the residual
 * shows that the factorization (or CG's tolerance) can do no better on this
 * system, and further steps would mostly measure rounding. Smaller stops
 * sooner; larger takes more steps that gain little, each costing a solve and
 * a multiply by K. A linear-solver rule, used only here. */
#define LIN_SYS_REFINE_CONTRACTION (0.5)

/** The solver's plain solve, in place: b = K^{-1} b. Direct solvers ignore s
 * and tol; CG solvers start from s (its first n entries; SCS_NULL: zero) and
 * stop at tol. Returns 0 on success. */
typedef scs_int (*ScsLinSysPlainSolve)(ScsLinSysWork *w, scs_float *b,
                                       const scs_float *s, scs_float tol);

/** Optional, for a solver with its own copy of K: r = b - K x with that K (as
 * of the last scs_update_lin_sys_diag_r), *nr = ||r||_inf, and *nb =
 * ||b||_inf unless nb is SCS_NULL. Writes r into the host vector r (n + m
 * entries, aliasing neither b nor x), unless the solver also has a correct
 * op: then r stays where the solver keeps its vectors for that op, and the
 * host r is not written. Returns 0 on success; on a failure the shared
 * routine prints that the product failed and the solve fails. */
typedef scs_int (*ScsLinSysResidual)(ScsLinSysWork *w, const scs_float *b,
                                     const scs_float *x, scs_float *r,
                                     scs_float *nr, scs_float *nb);

/** Optional, only with a residual op: out = x + K^{-1} r, where x and r are
 * those of the latest residual call (the shared routine always calls it right
 * after that residual, with the same x; out aliases neither x nor r). For a
 * solver whose vectors live elsewhere (cuDSS keeps r and x on the GPU, so it
 * neither copies r to the host nor x back). Returns 0 on success, and prints
 * why it failed otherwise. */
typedef scs_int (*ScsLinSysCorrect)(ScsLinSysWork *w, const scs_float *x,
                                    const scs_float *r, scs_float *out);

/** State of SCS(lin_sys_shared_solve), one per linear-solver workspace. */
typedef struct {
  /* K for SCS's products: the arrays given to scs_init_lin_sys_work, not
   * owned. SCS writes each new R into this same diag_r. */
  const ScsMatrix *A, *P;
  const scs_float *diag_r;
  scs_int n, m;
  /* 1: a CG solver (no comparison of the warm start with zero; refinement
   * restarts CG from x), 0: a direct solver */
  scs_int cg;
  ScsLinSysPlainSolve solve;
  ScsLinSysResidual residual; /* SCS_NULL: SCS's products */
  ScsLinSysCorrect correct;   /* SCS_NULL: x + the plain solve of a copy of r */
  /* n + m each: b0 a copy of the right-hand side (for residuals after the
   * solve has overwritten b; filled only when a residual will be needed), r
   * a residual, c a refinement candidate */
  scs_float *b0, *r, *c;
} ScsLinSysShared;

/**
 * Sets up sh for a solver's workspace w: call it at the end of
 * scs_init_lin_sys_work, with that function's A, P and diag_r, cg = 1 for a
 * CG solver (else 0), the solver's plain solve, and its residual and correct
 * ops or SCS_NULL. Allocates three vectors of n + m. Returns 0, or -1 if out
 * of memory (with a message; then call SCS(lin_sys_shared_free)).
 */
scs_int SCS(lin_sys_shared_init)(ScsLinSysShared *sh, const ScsMatrix *A,
                                 const ScsMatrix *P, const scs_float *diag_r,
                                 scs_int cg, ScsLinSysPlainSolve solve,
                                 ScsLinSysResidual residual,
                                 ScsLinSysCorrect correct);

/** Frees what SCS(lin_sys_shared_init) allocated (safe on a zeroed sh). */
void SCS(lin_sys_shared_free)(ScsLinSysShared *sh);

/**
 * scs_solve_lin_sys (linsys.h) for the solver with workspace w, whose sh was
 * set up by SCS(lin_sys_shared_init); the other arguments are those of
 * scs_solve_lin_sys, with its semantics.
 */
scs_int SCS(lin_sys_shared_solve)(ScsLinSysShared *sh, ScsLinSysWork *w,
                                  scs_float *b, const scs_float *s,
                                  scs_float tol, scs_int max_refine_steps,
                                  scs_float *res);

#ifdef __cplusplus
}
#endif
#endif
