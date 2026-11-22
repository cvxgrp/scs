/*
 * eigCG deflation, shared by the indirect (conjugate gradient) solvers: the
 * CPU one (linsys/cpu/indirect) and the GPU one (linsys/gpu/indirect).
 *
 * eigCG (Stathopoulos & Orginos 2010) observes the CG recurrence of the first
 * cold solve after each change of R (the solve for g = K^{-1} h, to
 * CG_BEST_TOL), without perturbing it, and extracts approximate eigenvectors
 * W for the smallest eigenvalues of the preconditioned reduced operator
 * G = R_x + P + A' R_y^{-1} A by thick-restarted Lanczos on CG's own scalars.
 * Every later solve, warm or cold, Galerkin-projects its initial residual
 * against span(W), which removes the slowly converging components; W spans
 * a near-invariant subspace, so CG does not reintroduce them. The set is
 * discarded at every scs_update_lin_sys_diag_r, and the next cold solve
 * harvests a new one.
 *
 * The algorithm (when to harvest, the window, the thick restart, the filter,
 * the Gram factor, the projection and its undo) is this one routine for both
 * solvers. Its small dense work (eigenproblems of the window's size by
 * LAPACK's syev, Gram-Schmidt, a Cholesky factorization) runs on the host,
 * and its work on vectors of length n through the solver's ScsEigcgOps,
 * where the solver keeps its vectors (on the GPU for the GPU solver). It
 * needs LAPACK and double precision: without USE_LAPACK, under SFLOAT, or
 * with DEFLATE_VECTORS 0 (glbopts.h), SCS(eigcg_init) leaves it off and
 * every other call does nothing.
 *
 * Memory. SCS(eigcg_init) allocates everything deflation uses, when the
 * solver is set up, and SCS(eigcg_free) frees it: W and AW (2 DEFLATE_VECTORS
 * vectors of length n) and the harvest's window (V, one vector per window
 * entry, and VS, 2 DEFLATE_VECTORS more), (4 DEFLATE_VECTORS + window) n
 * numbers in all, 210 n with the defaults, plus small arrays on the host. If
 * that cannot be allocated, SCS(eigcg_init) fails with a message and the
 * solver's setup fails with it. So whether a solve deflates is fixed when SCS
 * is built, and never depends on the memory free while it runs (on a shared
 * GPU, that depends on the other processes). The window is reused by every
 * harvest; a harvest reads only the window entries it has written itself
 * (SCS(eigcg_start) clears T, and each column of V and VS, and each small
 * array, is written before it is read), so its values are those it would
 * compute in fresh memory. (Before SCS 3.4.0 the CPU solver allocated the
 * window at each harvest and freed it after: the same values and the same
 * peak memory, held for less of the run, and plain CG when it did not fit.)
 *
 * Failures while solving. A failed syev or Cholesky factorization, or a
 * filter that keeps no vector, leaves that harvest without a set (plain CG
 * until the next harvest): deterministic functions of the values, as before
 * SCS 3.4.0. A failed operation of the solver (only the GPU solver's can
 * fail) fails the solve.
 *
 * A solver keeps one ScsEigcg in its workspace and calls:
 *  - SCS(eigcg_init) in scs_init_lin_sys_work (failing with it),
 *    SCS(eigcg_free) in scs_free_lin_sys_work, SCS(eigcg_reset) in
 *    scs_update_lin_sys_diag_r;
 *  - per plain solve, SCS(eigcg_begin_solve) first, and after CG
 *    SCS(eigcg_end_solve), whose return (the products of the extraction) it
 *    adds to its count of products; a solve may return between the two
 *    (a zero right-hand side), since the next begin starts afresh;
 *  - in CG, SCS(eigcg_deflate) on the initial residual (before the
 *    convergence test), SCS(eigcg_start) with the first preconditioned
 *    residual z and z'r, and SCS(eigcg_step) once per iteration, when
 *    alpha, the next beta, z and z'r are known and p is not yet updated.
 * This header is not installed (as linsys/scs_matrix.h).
 */

#ifndef EIGCG_H_GUARD
#define EIGCG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"

/* The safeguards of eigCG. They were set when deflation was written (August
 * 2026) without a sweep of their values.
 *
 * EIGCG_RITZ_RES_TOL: after a harvest, a Ritz pair (theta, w) is kept only if
 * its eigenresidual ||G w - theta M^{-1} w||_2 (M the Jacobi preconditioner)
 * is at most this fraction of ||theta M^{-1} w||_2. The projection's
 * coefficients scale like 1 / lambda, so the error of an inexact vector is
 * amplified by the smallest eigenvalue. Lower keeps fewer, more accurate
 * vectors; higher keeps rougher ones. Exact vectors pass at rounding level.
 *
 * EIGCG_UNDO_GROWTH: a projection that leaves the initial residual's norm
 * (CG_NORM) above this factor times its norm before is undone. For exact
 * eigenvectors the projection cannot grow the residual; growth means the
 * vectors' error is being amplified on this right-hand side. Lower undoes
 * more often, losing deflation's gain; higher lets an inexact set inject
 * error into the solve.
 *
 * EIGCG_MGS_RANK_TOL: in a thick restart, Gram-Schmidt (two passes) drops a
 * column of the restart stack whose norm is then below this: the numerical
 * rank of the stack, whose columns have unit norm. Lower keeps nearly
 * dependent columns, which the following eigenproblem then has to resolve;
 * higher drops columns that carry information.
 *
 * EIGCG_CHOL_MIN_PIVOT: the Cholesky factorization of W' G W fails, and the
 * harvest keeps no set, when a pivot is below this or not a number. It only
 * guards a division by zero: the filter above has already removed vectors
 * that would make W' G W nearly singular. */
#define EIGCG_RITZ_RES_TOL (0.25)
#define EIGCG_UNDO_GROWTH (1.5)
#define EIGCG_MGS_RANK_TOL (1e-8)
#define EIGCG_CHOL_MIN_PIVOT (1e-300)

/** The solver's operations on vectors of length n, for its ScsEigcg. Every
 * vector argument (named in upper case, or x, r, z, v, Kv) is storage from
 * alloc: n numbers, or n x cols stored by column (leading dimension n), where
 * the solver keeps its vectors (on the host for the CPU solver, on the GPU
 * for the GPU solver). Small matrices and coefficients (B, G, c) and scalars
 * are on the host. The ops that return scs_int return 0 on success; the CPU
 * solver's always succeed. Sizes are computed in size_t (n x cols can exceed
 * the range of scs_int). */
typedef struct {
  /* zeroed storage for rows x cols numbers, SCS_NULL if out of memory */
  scs_float *(*alloc)(ScsLinSysWork *w, scs_int rows, scs_int cols);
  /* frees storage from alloc (SCS_NULL is a no-op) */
  void (*free)(ScsLinSysWork *w, scs_float *v);
  /* dst = a * src */
  scs_int (*scale_copy)(ScsLinSysWork *w, scs_float *dst, const scs_float *src,
                        scs_float a);
  /* dst = src, cols columns of n (dst and src do not overlap) */
  scs_int (*copy)(ScsLinSysWork *w, scs_float *dst, const scs_float *src,
                  scs_int cols);
  /* C = A B: A is n x k, B is k x cols on the host (leading dimension k), C
   * is n x cols */
  scs_int (*mul_small)(ScsLinSysWork *w, scs_float *C, const scs_float *A,
                       scs_int k, const scs_float *B, scs_int cols);
  /* Kv = G v, then *rr = ||Kv - theta D v||_2^2 and *mm = ||theta D v||_2^2,
   * D = diag(G) = M^{-1} (the filter of EIGCG_RITZ_RES_TOL). The solver may
   * use its CG vectors as scratch: it runs after CG */
  scs_int (*probe)(ScsLinSysWork *w, const scs_float *v, scs_float *Kv,
                   scs_float theta, scs_float *rr, scs_float *mm);
  /* G = W' AW, k x k on the host: G[b * k + a] = W(:, a)' AW(:, b) */
  scs_int (*gram)(ScsLinSysWork *w, const scs_float *W, const scs_float *AW,
                  scs_int k, scs_float *G);
  /* c = W' r, k numbers on the host */
  scs_int (*coeffs)(ScsLinSysWork *w, const scs_float *W, scs_int k,
                    const scs_float *r, scs_float *c);
  /* x += W c and r -= AW c, or with undo nonzero x -= W c and r += AW c */
  scs_int (*project)(ScsLinSysWork *w, const scs_float *W, const scs_float *AW,
                     scs_int k, const scs_float *c, scs_float *x, scs_float *r,
                     scs_int undo);
  /* CG's convergence norm of r (CG_NORM) */
  scs_float (*norm)(ScsLinSysWork *w, const scs_float *r);
} ScsEigcgOps;

/** State of eigCG deflation, one per indirect solver's workspace. The
 * solver reads max, count, win (to size buffers of its own) and nothing
 * else. */
typedef struct {
  ScsLinSysWork *w;
  const ScsEigcgOps *ops;
  scs_int n;
  /* max: the most vectors kept (0: deflation off), count: the vectors of the
   * current set (0: none), harvest: this solve harvests */
  scs_int max, count, harvest;
  /* the set, from alloc: W (n x max), AW = G W (n x max); on the host: the
   * Cholesky factor of W' AW (max x max) and the projection's coefficients */
  scs_float *W, *AW, *G, *c;
  /* the window, reused by every harvest: from alloc, V the Lanczos (then
   * locked Ritz) basis (n x win) and VS the restart's scratch (n x 2 max);
   * on the host, T the projected operator (win x win), E the eigenvectors
   * (win x win), ev the eigenvalues (win), S the restart stack, TS its
   * scratch and Gr the restart transform (win x 2 max each), H the small
   * projection (2 max x 2 max), work syev's workspace (lwork) */
  scs_float *V, *VS, *T, *E, *ev, *S, *TS, *H, *Gr, *work;
  scs_int win, lwork;
  /* the harvest: j the vectors in the window, jd those whose row of T is
   * complete, step the CG iterations observed, have_prev whether pa and pb
   * (the previous alpha and beta) are set, dead whether this harvest stopped
   * observing (its window is still extracted), failed whether one of the
   * solver's ops failed (nothing is extracted, and the solve fails) */
  scs_int j, jd, step, have_prev, dead, failed;
  scs_float pa, pb;
} ScsEigcg;

/**
 * Sets up e for a solver with workspace w, n the length of its CG vectors,
 * and its ops (kept by pointer: give a static table). Deflation is on
 * (e->max > 0) with USE_LAPACK, without SFLOAT and with DEFLATE_VECTORS > 0;
 * then every array deflation uses is allocated here. Returns 0, or -1 (with
 * a message, and nothing left allocated) if that memory cannot be allocated
 * or LAPACK's workspace query fails; the solver's setup then fails.
 */
scs_int SCS(eigcg_init)(ScsEigcg *e, ScsLinSysWork *w, const ScsEigcgOps *ops,
                        scs_int n);

/** Frees all of e's memory and leaves deflation off (safe on a zeroed e,
 * and after a failed SCS(eigcg_init)). */
void SCS(eigcg_free)(ScsEigcg *e);

/** Discards the set (R changed: it lives in the old preconditioned
 * geometry); the next cold solve harvests a new one. */
void SCS(eigcg_reset)(ScsEigcg *e);

/** At the start of a plain solve: cold (no warm start) solves harvest when
 * no set is held. */
void SCS(eigcg_begin_solve)(ScsEigcg *e, scs_int cold);

/**
 * Projects CG's initial residual r against the set and adds the correction
 * to x (both vectors of the solver), undoing it if the residual grew by more
 * than EIGCG_UNDO_GROWTH. Costs no products. Does nothing without a set.
 * Returns 0, or -1 if one of the solver's ops failed (x and r are then
 * undefined).
 */
scs_int SCS(eigcg_deflate)(ScsEigcg *e, scs_float *x, scs_float *r);

/** CG's first preconditioned residual z and z'r (does nothing unless this
 * solve harvests). */
void SCS(eigcg_start)(ScsEigcg *e, const scs_float *z, scs_float ztr);

/** One CG iteration: alpha, the next beta, and the new z and z'r (does
 * nothing unless this solve harvests). */
void SCS(eigcg_step)(ScsEigcg *e, scs_float alpha, scs_float beta_next,
                     const scs_float *z, scs_float ztr);

/**
 * At the end of a plain solve: if it harvested, extracts the set from the
 * window. Returns the products spent (one multiply by G per vector probed),
 * 0 if it did not harvest, or -1 if one of the solver's ops failed during
 * this harvest (no set is kept; the solver fails the solve).
 */
scs_int SCS(eigcg_end_solve)(ScsEigcg *e);

#ifdef __cplusplus
}
#endif
#endif
