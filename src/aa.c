/*
 * Anderson acceleration.
 *
 * x: input iterate
 * x_prev: previous input iterate
 * f: f(x) output of map f applied to x
 * g: x - f (error)
 * g_prev: previous error
 * s: x - x_prev
 * y: g - g_prev
 * d: s - y = f - f_prev
 *
 * capital letters are the variables stacked columnwise
 * idx tracks current index where latest quantities written
 * idx cycles from left to right columns in matrix
 *
 * Type-I:
 * return f = f - (S - Y) * ( S'Y + r I)^{-1} ( S'g )
 *
 * Type-II:
 * return f = f - (S - Y) * ( Y'Y + r I)^{-1} ( Y'g )
 *
 * Both types reduce to the same regularized least-squares augmentation
 *     (A'B + r I) γ = A' g
 *     ⇔   [A; √r I]' [B; √r I] γ = [A; √r I]' [g; 0],
 * where A = S (type-I) or A = Y (type-II), and B = Y. We solve via a thin
 * QR factorization of the augmented A, which keeps the conditioning at
 * κ(A_aug) rather than the κ(A_aug)² that a normal-equations solve would
 * incur — critical near the optimum where Y rows are tiny and the Gram
 * matrix becomes numerically singular.
 */

#include <float.h>
#include <math.h>
#include <string.h>

#include "aa.h"
#include "scs_blas.h"

#ifndef SFLOAT
#define AA_EPS DBL_EPSILON
#else
#define AA_EPS FLT_EPSILON
#endif

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef USE_LAPACK

typedef void *ACCEL_WORK;

AaWork *aa_init(aa_int dim, aa_int mem, aa_int min_len, aa_int type1,
                aa_float regularization, aa_float relaxation,
                aa_float safeguard_factor, aa_float max_weight_norm,
                aa_int ir_max_steps, aa_int verbosity) {
  return SCS_NULL;
}
aa_float aa_apply(aa_float *f, const aa_float *x, AaWork *a) {
  return 0;
}
aa_int aa_safeguard(aa_float *f_new, aa_float *x_new, AaWork *a) {
  return 0;
}
void aa_finish(AaWork *a) {
}
void aa_reset(AaWork *a) {
}
AaStats aa_get_stats(const AaWork *a) {
  AaStats s;
  memset(&s, 0, sizeof(AaStats));
  s.last_aa_norm = NAN;
  return s;
}

#else

#if PROFILING > 0

#define TIME_TIC                                                               \
  timer __t;                                                                   \
  tic(&__t);
#define TIME_TOC toc(__func__, &__t);

#include <time.h>
typedef struct timer {
  struct timespec tic;
  struct timespec toc;
} timer;

void tic(timer *t) {
  clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

aa_float tocq(timer *t) {
  struct timespec temp;

  clock_gettime(CLOCK_MONOTONIC, &t->toc);

  if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
    temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
  } else {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
    temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
  }
  return (aa_float)temp.tv_sec * 1e3 + (aa_float)temp.tv_nsec / 1e6;
}

aa_float toc(const char *str, timer *t) {
  aa_float time = tocq(t);
  scs_printf("%s - time: %8.4f milli-seconds.\n", str, time);
  return time;
}

#else

#define TIME_TIC
#define TIME_TOC

#endif

#ifdef __cplusplus
extern "C" {
#endif

/* BLAS / LAPACK functions used */
aa_float BLAS(nrm2)(blas_int *n, aa_float *x, blas_int *incx);
void BLAS(axpy)(blas_int *n, aa_float *a, const aa_float *x, blas_int *incx,
                aa_float *y, blas_int *incy);
void BLAS(gemv)(const char *trans, const blas_int *m, const blas_int *n,
                const aa_float *alpha, const aa_float *a, const blas_int *lda,
                const aa_float *x, const blas_int *incx, const aa_float *beta,
                aa_float *y, const blas_int *incy);
void BLAS(gesv)(blas_int *n, blas_int *nrhs, aa_float *a, blas_int *lda,
                blas_int *ipiv, aa_float *b, blas_int *ldb, blas_int *info);
void BLAS(getrs)(const char *trans, const blas_int *n, const blas_int *nrhs,
                 const aa_float *a, const blas_int *lda, const blas_int *ipiv,
                 aa_float *b, const blas_int *ldb, blas_int *info);
void BLAS(scal)(const blas_int *n, const aa_float *a, aa_float *x,
                const blas_int *incx);
void BLAS(trsv)(const char *uplo, const char *trans, const char *diag,
                const blas_int *n, const aa_float *a, const blas_int *lda,
                aa_float *x, const blas_int *incx);
void BLAS(trmv)(const char *uplo, const char *trans, const char *diag,
                const blas_int *n, const aa_float *a, const blas_int *lda,
                aa_float *x, const blas_int *incx);
void BLAS(geqp3)(const blas_int *m, const blas_int *n, aa_float *a,
                 const blas_int *lda, blas_int *jpvt, aa_float *tau,
                 aa_float *work, const blas_int *lwork, blas_int *info);
void BLAS(ormqr)(const char *side, const char *trans, const blas_int *m,
                 const blas_int *n, const blas_int *k, const aa_float *a,
                 const blas_int *lda, const aa_float *tau, aa_float *c,
                 const blas_int *ldc, aa_float *work, const blas_int *lwork,
                 blas_int *info);

#ifdef __cplusplus
}
#endif

/* This file uses Anderson acceleration to improve the convergence of
 * a fixed point mapping.
 * At each iteration we solve a (small) regularized least-squares
 * problem via a pivoted QR factorization of an augmented matrix,
 * followed by iterative refinement on the reduced system.
 */

/* contains the necessary parameters to perform aa at each step */
struct ACCEL_WORK {
  aa_int type1;        /* bool, if true type 1 aa otherwise type 2 */
  aa_int mem;          /* aa memory */
  aa_int min_len;      /* min iterates before solve starts (1..mem) */
  aa_int dim;          /* variable dimension */
  aa_int iter;         /* current iteration */
  aa_int verbosity;    /* verbosity level, 0 is no printing */
  aa_int success;      /* was the last AA step successful or not */
  aa_int ir_max_steps; /* max iterative refinement passes, 0 disables */

  aa_float relaxation;       /* relaxation x and f, beta in some papers */
  aa_float regularization;   /* regularization */
  aa_float safeguard_factor; /* safeguard tolerance factor */
  aa_float max_weight_norm;  /* maximum norm of AA weights */

  aa_float *x;     /* x input to map*/
  aa_float *f;     /* f(x) output of map */
  aa_float *g;     /* x - f(x) */
  aa_float norm_g; /* ||x - f(x)|| */

  /* from previous iteration */
  aa_float *g_prev; /* x_prev - f(x_prev) */

  aa_float *Y; /* matrix of stacked y values */
  aa_float *S; /* matrix of stacked s values */
  aa_float *D; /* matrix of stacked d values = (S-Y) */

  /* Per-column cached L2 norms of S and Y (length mem). Each slot is
   * rewritten when its column is rewritten (update_accel_params), and
   * compute_regularization reduces these with a max-scaled sum to form
   * ||S||_F / ||Y||_F in O(mem) — versus the original O(dim·mem) nrm2
   * over the whole matrix. Storing per-column (not incremental "subtract
   * old / add new") avoids drift: near convergence Y columns can shrink
   * by many orders of magnitude and add/subtract updates fall below the
   * running sum's rounding floor, pegging the total above the true
   * Frobenius norm (observed on the κ=1e10 stress test). Storing the
   * norm (not its square) preserves nrm2's overflow/underflow safety —
   * squaring a 1e200 or 1e-200 column norm would hit ±inf / 0 even
   * though nrm2 on the whole matrix would produce a finite answer. */
  aa_float *nrm_s_col;
  aa_float *nrm_y_col;

  /* QR workspaces, sized for the augmented problem. */
  aa_float *A_aug;   /* (dim + mem) x mem  -- [A; √r I]; factored in place */
  aa_float *B_aug;   /* (dim + mem) x mem  -- [Y; √r I] (type-I only) */
  aa_float *c_aug;   /* (dim + mem)        -- [g; 0], overwritten by Q' c */
  aa_float *tau;     /* mem                -- Householder scalars */
  aa_float *qr_work; /* lwork              -- LAPACK scratch for geqp3/ormqr */
  blas_int qr_lwork; /* size of qr_work, chosen via workspace query at init */
  blas_int *jpvt;    /* mem                -- column permutation from geqp3 */

  aa_float *W;      /* mem x mem scratch: Q' B_aug top block (type-I gesv) */
  aa_float *W_orig; /* mem x mem: W copy preserved for iterative refinement */
  blas_int *ipiv;   /* gesv permutation (type-I) */

  /* Iterative refinement scratches (size mem each). Separate from `work`
   * so the IR dance doesn't clobber γ or the safeguard scratch. */
  aa_float *gamma_red;   /* permuted/truncated γ (length `rank`) */
  aa_float *c_top_save;  /* original RHS preserved across the solve */
  aa_float *ir_res;      /* residual/correction vector */

  /* Dual-use scratch: dim-sized buffer for aa_safeguard's x_new - f_new
   * diff, and also the natural-order γ home (len ≤ mem entries) inside
   * solve(). Allocated as max(mem, dim) so both uses fit. */
  aa_float *work;

  aa_float *x_work; /* workspace (= x) for when relaxation != 1.0 */

  /* Lifetime diagnostics (see AaStats in aa_stats.h). NOT cleared by
   * aa_reset — the internal reset path fires on safeguard rejection,
   * and you want the rejection to stay visible in the counters. */
  aa_int n_accept;
  aa_int n_reject_lapack;
  aa_int n_reject_rank0;
  aa_int n_reject_nonfinite;
  aa_int n_reject_weight_cap;
  aa_int n_safeguard_reject;
  aa_int last_rank;
  aa_float last_aa_norm;
  aa_float last_regularization;
};

/* Reduce a length-`mem` vector of nonnegative column norms to a single
 * Frobenius norm ||A||_F = sqrt(Σ nrm_col_i²), using the classic
 * max-scaled sum of squares so we don't square-then-overflow on a
 * large column or square-then-underflow on a tiny one. This mirrors
 * what nrm2 does internally across elements, but here across column
 * norms — each nrm_col_i is already nrm2-safe per column. */
static aa_float frob_from_col_norms(const aa_float *nrm_col, aa_int mem) {
  aa_int i;
  aa_float m = 0;
  for (i = 0; i < mem; ++i) {
    if (nrm_col[i] > m) m = nrm_col[i];
  }
  if (m == 0) return 0;
  aa_float sumsq = 0;
  for (i = 0; i < mem; ++i) {
    aa_float t = nrm_col[i] / m;
    sumsq += t * t;
  }
  return m * sqrt(sumsq);
}

/* Tikhonov regularization scaled with the problem. Matches the prior
 * behavior's intent (r grows with the magnitude of A'B so `regularization`
 * stays unitless), but uses the cheap Frobenius-norm upper bound
 *     ||A'B||_F ≤ ||A||_F · ||B||_F
 * instead of maintaining a Gram matrix. For type-II A == B so this is
 * ||Y||_F², the same scale as the previous ||Y'Y||_F up to a factor
 * ≤ √mem. Reduces the per-column cached norms (O(mem)) rather than a
 * fresh nrm2 over dim·mem entries. */
static aa_float compute_regularization(AaWork *a) {
  TIME_TIC
  aa_float nrm_y = frob_from_col_norms(a->nrm_y_col, a->mem);
  aa_float nrm_a = a->type1 ? frob_from_col_norms(a->nrm_s_col, a->mem) : nrm_y;
  aa_float r = a->regularization * nrm_a * nrm_y;
  if (a->verbosity > 2) {
    scs_printf("iter: %i, ||A||_F %.2e, ||Y||_F %.2e, r: %.2e\n",
               (int)a->iter, nrm_a, nrm_y, r);
  }
  TIME_TOC
  return r;
}

/* Build [M; √r I_len] column-major into `dst` with fixed leading dim
 * (dim + mem). When len < mem we zero-pad the unused trailing rows so
 * the QR factorization still operates on a well-defined (dim+mem) x len
 * block; the zero rows don't change the solve. */
static void build_augmented(aa_float *dst, const aa_float *src, aa_int dim,
                            aa_int mem, aa_int len, aa_float sqrt_r) {
  aa_int i;
  aa_int aug_rows = dim + mem;
  for (i = 0; i < len; ++i) {
    aa_float *col = &dst[i * aug_rows];
    memcpy(col, &src[i * dim], dim * sizeof(aa_float));
    memset(&col[dim], 0, mem * sizeof(aa_float));
    col[dim + i] = sqrt_r;
  }
}

/* initialize accel params, in particular x_prev, f_prev, g_prev */
static void init_accel_params(const aa_float *x, const aa_float *f, AaWork *a) {
  TIME_TIC
  blas_int bdim = (blas_int)a->dim;
  aa_float neg_onef = -1.0;
  blas_int one = 1;
  /* x_prev = x */
  memcpy(a->x, x, sizeof(aa_float) * a->dim);
  /* f_prev = f */
  memcpy(a->f, f, sizeof(aa_float) * a->dim);
  /* g_prev = x */
  memcpy(a->g_prev, x, sizeof(aa_float) * a->dim);
  /* g_prev = x_prev - f_prev */
  BLAS(axpy)(&bdim, &neg_onef, f, &one, a->g_prev, &one);
  TIME_TOC
}

/* updates the workspace parameters for aa for this iteration
 *
 * Writes this iteration's s, d, y columns directly into S, D, Y at slot
 * `idx` — no intermediate scratch. Numerically sensitive because:
 *
 *  - y is computed as g - g_prev (ONE rounding into a cancellation-prone
 *    quantity). Deriving y from s - d would add two extra roundings and
 *    make y noticeably worse near convergence where g and g_prev are
 *    tiny and nearly equal.
 *
 *  - The reads of a->x, a->f, a->g_prev all require the PREVIOUS
 *    iteration's values, so state advance (x_prev <- x, f_prev <- f,
 *    g_prev <- g) must happen AFTER everything that reads them. s uses
 *    old a->x; d uses old a->f; y uses old a->g_prev. */
static void update_accel_params(const aa_float *x, const aa_float *f, AaWork *a,
                                aa_int len) {
  /* Entry invariant:  a->x == x_prev, a->f == f_prev, a->g_prev == g_prev. */
  TIME_TIC
  aa_int idx = (a->iter - 1) % a->mem;
  blas_int one = 1;
  blas_int bdim = (blas_int)a->dim;
  aa_float neg_onef = -1.0;
  aa_float *s_col = &(a->S[idx * a->dim]);
  aa_float *d_col = &(a->D[idx * a->dim]);
  aa_float *y_col = &(a->Y[idx * a->dim]);

  /* S[:, idx] = x - x_prev  (reads old a->x). */
  memcpy(s_col, x, sizeof(aa_float) * a->dim);
  BLAS(axpy)(&bdim, &neg_onef, a->x, &one, s_col, &one);

  /* D[:, idx] = f - f_prev  (reads old a->f). */
  memcpy(d_col, f, sizeof(aa_float) * a->dim);
  BLAS(axpy)(&bdim, &neg_onef, a->f, &one, d_col, &one);

  /* g = x - f  (this iteration's residual; needed for the solve RHS). */
  memcpy(a->g, x, sizeof(aa_float) * a->dim);
  BLAS(axpy)(&bdim, &neg_onef, f, &one, a->g, &one);

  /* Y[:, idx] = g - g_prev  (reads old a->g_prev; single-rounding y). */
  memcpy(y_col, a->g, sizeof(aa_float) * a->dim);
  BLAS(axpy)(&bdim, &neg_onef, a->g_prev, &one, y_col, &one);

  /* Update the per-column cached norms for the slot we just rewrote.
   * compute_regularization reduces these on demand. Store the norm
   * itself (not its square) — squaring here would throw away the
   * overflow/underflow safety that nrm2 guarantees. */
  a->nrm_s_col[idx] = BLAS(nrm2)(&bdim, s_col, &one);
  a->nrm_y_col[idx] = BLAS(nrm2)(&bdim, y_col, &one);

  /* State advance for next iter: (x_prev, f_prev, g_prev) <- (x, f, g).
   * Must follow all the reads above. */
  memcpy(a->x, x, sizeof(aa_float) * a->dim);
  memcpy(a->f, f, sizeof(aa_float) * a->dim);
  memcpy(a->g_prev, a->g, sizeof(aa_float) * a->dim);

  /* Relaxation scratch (mirror of x); only present when relaxation != 1.0. */
  if (a->x_work) {
    memcpy(a->x_work, x, sizeof(aa_float) * a->dim);
  }

  /* ||g|| = ||x - f|| (current residual norm, used by the safeguard). */
  a->norm_g = BLAS(nrm2)(&bdim, a->g, &one);

  TIME_TOC
}

/* f = (1-relaxation) * \sum_i a_i x_i + relaxation * \sum_i a_i f_i */
static void relax(aa_float *f, AaWork *a, aa_int len) {
  TIME_TIC
  /* x_work = x initially */
  blas_int bdim = (blas_int)(a->dim), one = 1, blen = (blas_int)len;
  aa_float onef = 1.0, neg_onef = -1.0;
  aa_float one_m_relaxation = 1. - a->relaxation;
  /* x_work = x - S * work */
  BLAS(gemv)
  ("NoTrans", &bdim, &blen, &neg_onef, a->S, &bdim, a->work, &one, &onef,
   a->x_work, &one);
  /* f = relaxation * f */
  BLAS(scal)(&bdim, &a->relaxation, f, &one);
  /* f += (1 - relaxation) * x_work */
  BLAS(axpy)(&bdim, &one_m_relaxation, a->x_work, &one, f, &one);
  TIME_TOC
}

/* Solve the regularized normal equations (A'B + rI) γ = A'g via a
 * pivoted QR (geqp3) of the augmented matrix [A; √r I]. Column pivoting
 * exposes the numerical rank directly in the diagonal of R: we truncate
 * at the first diagonal whose magnitude falls below len·ε·|R_11|
 * (dim-independent: the inner LS is a len-column problem, so its noise
 * floor scales with column count, not the caller's state dim) and solve
 * the smaller, well-conditioned system (graceful degradation instead of
 * hard reset on near-rank-deficiency). Iterative refinement on the
 * reduced system recovers digits lost to gesv/trsv rounding; the loop
 * auto-stops when the correction no longer contracts and is capped at
 * ir_max_steps (see aa_init). γ is then validated against
 * max_weight_norm in the L2 sense. */
static aa_float solve(aa_float *f, AaWork *a, aa_int len) {
  TIME_TIC
  blas_int info = -1, bdim = (blas_int)(a->dim), one = 1, blen = (blas_int)len;
  /* Leading dim is fixed to dim+mem regardless of len so the buffers
   * allocated in aa_init match the strides used here. Unused trailing
   * rows are kept zero by build_augmented. */
  blas_int aug_rows = bdim + (blas_int)a->mem;
  blas_int bmem = (blas_int)a->mem;
  aa_float onef = 1.0, neg_onef = -1.0, aa_norm;
  aa_float *A_src = a->type1 ? a->S : a->Y;
  aa_float *gamma = a->work; /* natural-order γ, len entries used by gemv below */
  aa_int i;
  aa_int rank = 0;
  blas_int brank;

  /* Three regularization modes:
   *   regularization > 0  : problem-scaled   r = regularization * ||A||_F ||Y||_F
   *   regularization < 0  : pinned absolute  r = -regularization     (Frobenius skipped)
   *   regularization == 0 : unregularized    r = 0
   * Pinned mode gives a knob for applications where the problem scale is
   * known and the caller wants a stable, scale-free regularizer. */
  aa_float r;
  if (a->regularization > 0) {
    r = compute_regularization(a);
  } else if (a->regularization < 0) {
    r = -a->regularization;
  } else {
    r = 0.0;
  }
  aa_float sqrt_r = (r > 0) ? sqrt(r) : 0.0;

  /* 1. Build A_aug = [A; √r I_len]; factor with column pivoting. geqp3
   *    requires jpvt zeroed on entry so it is free to choose the pivot
   *    order (nonzero entries would be treated as user-forced pivots). */
  build_augmented(a->A_aug, A_src, a->dim, a->mem, len, sqrt_r);
  for (i = 0; i < len; ++i) a->jpvt[i] = 0;
  BLAS(geqp3)(&aug_rows, &blen, a->A_aug, &aug_rows, a->jpvt, a->tau,
              a->qr_work, &a->qr_lwork, &info);
  /* Capture geqp3's info before the rank-0 path below overwrites it; the
   * reject-cause attribution needs to distinguish a genuine LAPACK failure
   * from "the matrix went numerically to zero." */
  blas_int lapack_info = info;

  /* 2. Rank estimation. geqp3 sorts |R_ii| non-increasingly; find the
   *    largest `rank` with |R_rank-1,rank-1| ≥ tol. A rank of zero means
   *    the whole Ã is numerically zero — hand off to aa_reset below. */
  if (info == 0) {
    aa_float r11 = fabs(a->A_aug[0]);
    if (r11 > 0) {
      /* Column-count-based rank tolerance: the effective LS problem has
       * `len` columns, so the rounding floor for rank determination
       * scales with `len`, not (dim + mem). Decoupling from `dim` avoids
       * falsely dropping columns that are healthy relative to the
       * regularizer at large state dimensions. */
      aa_float tol = r11 * (aa_float)len * AA_EPS;
      for (rank = 0; rank < len; ++rank) {
        if (fabs(a->A_aug[rank * aug_rows + rank]) < tol) break;
      }
    }
    if (rank == 0) info = 1;
  }
  brank = (blas_int)rank;

  /* 3. c_aug = [g; 0]; overwrite with Q' c_aug. We only need the first
   *    `rank` entries (= Q_rank' c̃); pass `rank` to ormqr so it applies
   *    only the reflectors we care about. */
  if (info == 0) {
    memcpy(a->c_aug, a->g, a->dim * sizeof(aa_float));
    memset(&a->c_aug[a->dim], 0, a->mem * sizeof(aa_float));
    BLAS(ormqr)
    ("Left", "Trans", &aug_rows, &one, &brank, a->A_aug, &aug_rows, a->tau,
     a->c_aug, &aug_rows, a->qr_work, &a->qr_lwork, &info);
  }

  /* 4. Solve the reduced rank×rank system for γ_red (pivoted order),
   *    then un-permute into the natural-order γ consumed by `f -= D γ`. */
  if (info == 0) {
    /* Preserve the RHS for iterative refinement below. */
    memcpy(a->c_top_save, a->c_aug, rank * sizeof(aa_float));

    if (a->type1) {
      /* Type-I: build B_aug with the pivoted Y columns (first `rank` only),
       * apply Q', extract top-left rank×rank block into W, solve
       * W γ_red = c_top. The √rI block is reshuffled too — column i of
       * the permuted B̃ has √r at row (jpvt[i]-1). */
      for (i = 0; i < rank; ++i) {
        aa_int piv = a->jpvt[i] - 1; /* LAPACK jpvt is 1-indexed */
        aa_float *col = &a->B_aug[i * aug_rows];
        memcpy(col, &a->Y[piv * a->dim], a->dim * sizeof(aa_float));
        memset(&col[a->dim], 0, a->mem * sizeof(aa_float));
        col[a->dim + piv] = sqrt_r;
      }
      BLAS(ormqr)
      ("Left", "Trans", &aug_rows, &brank, &brank, a->A_aug, &aug_rows,
       a->tau, a->B_aug, &aug_rows, a->qr_work, &a->qr_lwork, &info);
      if (info == 0) {
        /* W (mem×mem, LDA=mem) holds the rank×rank top-left block. */
        for (i = 0; i < rank; ++i) {
          memcpy(&a->W[i * a->mem], &a->B_aug[i * aug_rows],
                 rank * sizeof(aa_float));
          memcpy(&a->W_orig[i * a->mem], &a->W[i * a->mem],
                 rank * sizeof(aa_float));
        }
        memcpy(a->gamma_red, a->c_top_save, rank * sizeof(aa_float));
        BLAS(gesv)
        (&brank, &one, a->W, &bmem, a->ipiv, a->gamma_red, &brank, &info);
        if (info == 0) {
          /* Iterative refinement: repeat while δ is still contracting,
           * capped at ir_max_steps. Each step: ρ = c_top - W_orig γ_red,
           * solve W δ = ρ with the LU already in W, γ_red += δ. Stop when
           * ‖δ_k‖ ≥ 0.5·‖δ_{k-1}‖ (we've hit the working-precision floor
           * and further steps won't help). */
          aa_float prev_dnorm = 0.0;
          aa_int k;
          for (k = 0; k < a->ir_max_steps; ++k) {
            aa_float dnorm;
            memcpy(a->ir_res, a->c_top_save, rank * sizeof(aa_float));
            BLAS(gemv)
            ("NoTrans", &brank, &brank, &neg_onef, a->W_orig, &bmem,
             a->gamma_red, &one, &onef, a->ir_res, &one);
            BLAS(getrs)
            ("NoTrans", &brank, &one, a->W, &bmem, a->ipiv, a->ir_res,
             &brank, &info);
            if (info != 0) break;
            dnorm = BLAS(nrm2)(&brank, a->ir_res, &one);
            BLAS(axpy)(&brank, &onef, a->ir_res, &one, a->gamma_red, &one);
            if (k > 0 && dnorm >= 0.5 * prev_dnorm) break;
            prev_dnorm = dnorm;
          }
        }
      }
    } else {
      /* Type-II: B̃ = Ã, so Q' B̃ = R. Solve R u = c_top for the permuted
       * solution u; the rank×rank leading block of R lives in the upper
       * triangle of A_aug. Rank truncation above already guaranteed the
       * diagonal is nonzero through index rank-1. */
      memcpy(a->gamma_red, a->c_top_save, rank * sizeof(aa_float));
      BLAS(trsv)("Upper", "NoTrans", "NonUnit", &brank, a->A_aug,
                 &aug_rows, a->gamma_red, &one);
      /* Iterative refinement, capped at ir_max_steps with early stop when
       * δ stops contracting. ρ = c_top - R u via trmv + subtract; solve
       * R δ = ρ via trsv; u += δ. */
      {
        aa_float prev_dnorm = 0.0;
        aa_int k;
        for (k = 0; k < a->ir_max_steps; ++k) {
          aa_float dnorm;
          memcpy(a->ir_res, a->gamma_red, rank * sizeof(aa_float));
          BLAS(trmv)("Upper", "NoTrans", "NonUnit", &brank, a->A_aug,
                     &aug_rows, a->ir_res, &one);
          for (i = 0; i < rank; ++i) {
            a->ir_res[i] = a->c_top_save[i] - a->ir_res[i];
          }
          BLAS(trsv)("Upper", "NoTrans", "NonUnit", &brank, a->A_aug,
                     &aug_rows, a->ir_res, &one);
          dnorm = BLAS(nrm2)(&brank, a->ir_res, &one);
          BLAS(axpy)(&brank, &onef, a->ir_res, &one, a->gamma_red, &one);
          if (k > 0 && dnorm >= 0.5 * prev_dnorm) break;
          prev_dnorm = dnorm;
        }
      }
    }

    /* Un-permute γ_red into γ (natural column order); zero the rest so
     * the `f -= D γ` gemv below sees a well-defined full-length vector. */
    if (info == 0) {
      memset(gamma, 0, len * sizeof(aa_float));
      for (i = 0; i < rank; ++i) {
        gamma[a->jpvt[i] - 1] = a->gamma_red[i];
      }
    }
  }

  /* 5. Validate γ via ‖γ‖₂ against max_weight_norm. */
  aa_norm = (info == 0) ? BLAS(nrm2)(&blen, gamma, &one) : -1.0;

  /* Record diagnostics for this solve, regardless of accept/reject.
   * NaN last_aa_norm signals "no valid norm this solve" — distinguishing
   * the genuine-zero case (rank collapse gives aa_norm = 0 legitimately)
   * from a failed/rejected solve. */
  a->last_rank = rank;
  a->last_regularization = r;
  a->last_aa_norm = (info == 0 && isfinite(aa_norm)) ? aa_norm : NAN;

  if (a->verbosity > 1) {
    scs_printf("AA type %i, iter: %i, len %i, rank %i, info: %i, aa_norm %.2e\n",
               a->type1 ? 1 : 2, (int)a->iter, (int)len, (int)rank, (int)info,
               aa_norm);
  }

  if (info != 0 || !isfinite(aa_norm) || aa_norm >= a->max_weight_norm) {
    if (a->verbosity > 0) {
      scs_printf("Error in AA type %i, iter: %i, len %i, rank %i, info: %i, "
                 "aa_norm %.2e\n",
                 a->type1 ? 1 : 2, (int)a->iter, (int)len, (int)rank, (int)info,
                 aa_norm);
    }
    /* Attribute the rejection to exactly one cause, in priority order.
     * lapack_info is the original geqp3 return; the rank-0 path above may
     * have set info=1 but that is the bookkeeping trick, not a LAPACK
     * failure. Without this ordering, rank-0 would be miscounted as
     * "lapack" via info. */
    if (lapack_info != 0) {
      a->n_reject_lapack++;
    } else if (rank == 0) {
      a->n_reject_rank0++;
    } else if (!isfinite(aa_norm)) {
      a->n_reject_nonfinite++;
    } else {
      a->n_reject_weight_cap++;
    }
    a->success = 0;
    aa_reset(a);
    TIME_TOC
    if (!isfinite(aa_norm)) aa_norm = -1.0;
    return (aa_norm < 0) ? aa_norm : -aa_norm;
  }

  /* f -= D γ */
  BLAS(gemv)
  ("NoTrans", &bdim, &blen, &neg_onef, a->D, &bdim, gamma, &one, &onef, f,
   &one);

  if (a->relaxation != 1.0) {
    relax(f, a, len);
  }

  a->success = 1;
  TIME_TOC
  return aa_norm;
}

/*
 * API functions below this line, see aa.h for descriptions.
 */
AaWork *aa_init(aa_int dim, aa_int mem, aa_int min_len, aa_int type1,
                aa_float regularization, aa_float relaxation,
                aa_float safeguard_factor, aa_float max_weight_norm,
                aa_int ir_max_steps, aa_int verbosity) {
  TIME_TIC
  AaWork *a;
  aa_int mem_clamped = MIN(mem, dim);
  /* `regularization` is accepted with either sign: positive = scaled by
   * ||A||_F ||Y||_F; negative = pinned absolute |regularization|; zero = off.
   * Only NaN / non-finite values are rejected (via the !isfinite check).
   * min_len < 1 is rejected when mem > 0; min_len > mem_clamped is
   * silently clamped down — same treatment the `mem` argument already
   * gets against `dim`, so callers can pass `min_len = mem` without
   * caring whether mem exceeded dim. When mem == 0 (AA off), min_len
   * is ignored entirely. */
  if (dim <= 0 || mem < 0 || !isfinite(regularization) ||
      relaxation < 0 || relaxation > 2 ||
      safeguard_factor < 0 || max_weight_norm <= 0 ||
      ir_max_steps < 0 ||
      (mem_clamped > 0 && min_len < 1)) {
    scs_printf("Invalid AA parameters.\n");
    return SCS_NULL;
  }
  a = (AaWork *)scs_calloc(1, sizeof(AaWork));
  if (!a) {
    scs_printf("Failed to allocate memory for AA.\n");
    return SCS_NULL;
  }
  a->type1 = type1;
  a->iter = 0;
  a->dim = dim;
  a->mem = mem_clamped; /* clamped to dim for rank stability */
  if (mem > dim && verbosity > 0) {
    scs_printf("AA: mem (%d) > dim (%d); clamping mem to dim.\n",
               (int)mem, (int)dim);
  }
  a->min_len = mem_clamped > 0 ? MIN(min_len, mem_clamped) : 0;
  a->regularization = regularization;
  a->relaxation = relaxation;
  a->safeguard_factor = safeguard_factor;
  a->max_weight_norm = max_weight_norm;
  a->ir_max_steps = ir_max_steps;
  a->success = 0;
  a->verbosity = verbosity;
  /* Counters are already zero from calloc; only last_aa_norm needs an
   * explicit sentinel so callers can distinguish "never solved" from a
   * legitimate zero norm (which never happens on a successful solve, but
   * 0 is a bad signal either way). */
  a->last_aa_norm = NAN;
  if (a->mem <= 0) {
    return a;
  }

  a->x = (aa_float *)scs_calloc(a->dim, sizeof(aa_float));
  a->f = (aa_float *)scs_calloc(a->dim, sizeof(aa_float));
  a->g = (aa_float *)scs_calloc(a->dim, sizeof(aa_float));

  a->g_prev = (aa_float *)scs_calloc(a->dim, sizeof(aa_float));

  a->Y = (aa_float *)scs_calloc(a->dim * a->mem, sizeof(aa_float));
  a->S = (aa_float *)scs_calloc(a->dim * a->mem, sizeof(aa_float));
  a->D = (aa_float *)scs_calloc(a->dim * a->mem, sizeof(aa_float));

  {
    aa_int aug_rows = a->dim + a->mem;
    a->A_aug = (aa_float *)scs_calloc((size_t)aug_rows * a->mem, sizeof(aa_float));
    a->c_aug = (aa_float *)scs_calloc((size_t)aug_rows, sizeof(aa_float));
    a->tau = (aa_float *)scs_calloc(a->mem, sizeof(aa_float));
    a->jpvt = (blas_int *)scs_calloc(a->mem, sizeof(blas_int));

    /* Scratches for iterative refinement (both types). */
    a->gamma_red = (aa_float *)scs_calloc(a->mem, sizeof(aa_float));
    a->c_top_save = (aa_float *)scs_calloc(a->mem, sizeof(aa_float));
    a->ir_res = (aa_float *)scs_calloc(a->mem, sizeof(aa_float));

    /* Per-column cached norms used by compute_regularization. */
    a->nrm_s_col = (aa_float *)scs_calloc(a->mem, sizeof(aa_float));
    a->nrm_y_col = (aa_float *)scs_calloc(a->mem, sizeof(aa_float));

    /* type-I needs a second augmented buffer and mem×mem gesv scratches;
     * W_orig preserves W across gesv so iterative refinement can form
     * the residual c_top − W γ. */
    if (type1) {
      a->B_aug = (aa_float *)scs_calloc((size_t)aug_rows * a->mem, sizeof(aa_float));
      a->W = (aa_float *)scs_calloc((size_t)a->mem * a->mem, sizeof(aa_float));
      a->W_orig = (aa_float *)scs_calloc((size_t)a->mem * a->mem, sizeof(aa_float));
      a->ipiv = (blas_int *)scs_calloc(a->mem, sizeof(blas_int));
    } else {
      a->B_aug = SCS_NULL;
      a->W = SCS_NULL;
      a->W_orig = SCS_NULL;
      a->ipiv = SCS_NULL;
    }

    a->work = (aa_float *)scs_calloc(MAX(a->mem, a->dim), sizeof(aa_float));
    if (relaxation != 1.0) {
      a->x_work = (aa_float *)scs_calloc(a->dim, sizeof(aa_float));
    } else {
      a->x_work = SCS_NULL;
    }

    /* Check every allocation before the LAPACK workspace query below. The
     * query passes A_aug/jpvt/tau/c_aug/B_aug into geqp3/ormqr; if any of
     * them is NULL we'd dereference inside LAPACK instead of returning
     * cleanly. qr_work is still NULL here — it's allocated after the query. */
    if (!a->x || !a->f || !a->g || !a->g_prev ||
        !a->Y || !a->S || !a->D ||
        !a->A_aug || !a->c_aug || !a->tau || !a->jpvt ||
        !a->gamma_red || !a->c_top_save || !a->ir_res ||
        !a->nrm_s_col || !a->nrm_y_col ||
        (type1 && (!a->B_aug || !a->W || !a->W_orig || !a->ipiv)) ||
        !a->work ||
        (relaxation != 1.0 && !a->x_work)) {
      scs_printf("Failed to allocate memory for AA.\n");
      aa_finish(a);
      return SCS_NULL;
    }

    /* LAPACK workspace query: ask geqp3 and ormqr for their preferred lwork,
     * then take the max. lwork = -1 makes the routine write the optimal
     * size into work[0] without doing any factoring. geqp3 typically wants
     * more scratch than geqrf because it also maintains column norms. The
     * optimal size is returned in an aa_float slot; round up with ceil
     * before casting so a value like 255.9999 doesn't truncate to 255 and
     * under-allocate. */
    {
      blas_int b_aug = (blas_int)aug_rows;
      blas_int b_mem = (blas_int)a->mem;
      blas_int b_neg_one = -1;
      blas_int info_q = 0;
      aa_float q_geqp3 = 0.0, q_ormqr_c = 0.0, q_ormqr_b = 0.0;
      BLAS(geqp3)(&b_aug, &b_mem, a->A_aug, &b_aug, a->jpvt, a->tau,
                  &q_geqp3, &b_neg_one, &info_q);
      {
        blas_int b_one = 1;
        BLAS(ormqr)
        ("Left", "Trans", &b_aug, &b_one, &b_mem, a->A_aug, &b_aug, a->tau,
         a->c_aug, &b_aug, &q_ormqr_c, &b_neg_one, &info_q);
      }
      if (type1) {
        BLAS(ormqr)
        ("Left", "Trans", &b_aug, &b_mem, &b_mem, a->A_aug, &b_aug, a->tau,
         a->B_aug, &b_aug, &q_ormqr_b, &b_neg_one, &info_q);
      }
      {
        aa_float lwork_f = q_geqp3;
        if (q_ormqr_c > lwork_f) lwork_f = q_ormqr_c;
        if (q_ormqr_b > lwork_f) lwork_f = q_ormqr_b;
        /* Floor at mem — some LAPACK builds return modest sizes; keep
         * a sane minimum. calloc of zero is implementation-defined. */
        if (lwork_f < (aa_float)a->mem) lwork_f = (aa_float)a->mem;
        a->qr_lwork = (blas_int)ceil(lwork_f);
        a->qr_work = (aa_float *)scs_calloc((size_t)a->qr_lwork, sizeof(aa_float));
      }
    }
    if (!a->qr_work) {
      scs_printf("Failed to allocate memory for AA.\n");
      aa_finish(a);
      return SCS_NULL;
    }
  }
  TIME_TOC
  return a;
}

aa_float aa_apply(aa_float *f, const aa_float *x, AaWork *a) {
  TIME_TIC
  aa_float aa_norm = 0;
  aa_int len = MIN(a->iter, a->mem);
  a->success = 0; /* if we make an AA step we set this to 1 later */
  if (a->mem <= 0) {
    TIME_TOC
    return aa_norm; /* 0 */
  }
  if (a->iter == 0) {
    /* if first iteration then seed params for next iter */
    init_accel_params(x, f, a);
    a->iter++;
    TIME_TOC
    return aa_norm; /* 0 */
  }
  /* set various accel quantities */
  update_accel_params(x, f, a, len);

  /* Hold off the solve until we have min_len residual pairs buffered. */
  if (a->iter >= a->min_len) {
    /* solve linear system, new point overwrites f if successful.
     * Rejection causes are counted inside solve() where the specific
     * failure mode is known; here we only count acceptances. */
    aa_norm = solve(f, a, len);
    if (aa_norm > 0) {
      a->n_accept++;
    }
  }
  a->iter++;
  TIME_TOC
  return aa_norm;
}

aa_int aa_safeguard(aa_float *f_new, aa_float *x_new, AaWork *a) {
  TIME_TIC
  blas_int bdim = (blas_int)a->dim;
  blas_int one = 1;
  aa_float neg_onef = -1.0;
  aa_float norm_diff;
  if (a->mem <= 0) {
    /* degenerate workspace, nothing to safeguard against */
    TIME_TOC
    return 0;
  }
  if (!a->success) {
    /* last AA update was not successful, no need for safeguarding */
    TIME_TOC
    return 0;
  }

  /* reset success indicator in case safeguarding called multiple times */
  a->success = 0;

  /* NB: a->work is used here as a dim-sized scratch, but elsewhere (in solve)
   * only as a len-sized (<=mem) scratch. This is why it is allocated with
   * MAX(mem, dim) in aa_init — do not shrink it to mem. */
  /* work = x_new */
  memcpy(a->work, x_new, a->dim * sizeof(aa_float));
  /* work = x_new - f_new */
  BLAS(axpy)(&bdim, &neg_onef, f_new, &one, a->work, &one);
  /* norm_diff = || f_new - x_new || */
  norm_diff = BLAS(nrm2)(&bdim, a->work, &one);
  /* g = f - x */
  if (norm_diff > a->safeguard_factor * a->norm_g) {
    /* in this case we reject the AA step and reset */
    memcpy(f_new, a->f, a->dim * sizeof(aa_float));
    memcpy(x_new, a->x, a->dim * sizeof(aa_float));
    if (a->verbosity > 0) {
      scs_printf("AA rejection, iter: %i, norm_diff %.4e, prev_norm_diff %.4e\n",
                 (int)a->iter, norm_diff, a->norm_g);
    }
    a->n_safeguard_reject++;
    aa_reset(a);
    TIME_TOC
    return -1;
  }
  TIME_TOC
  return 0;
}

void aa_finish(AaWork *a) {
  if (a) {
    scs_free(a->x);
    scs_free(a->f);
    scs_free(a->g);
    scs_free(a->g_prev);
    scs_free(a->Y);
    scs_free(a->S);
    scs_free(a->D);
    scs_free(a->A_aug);
    scs_free(a->B_aug);
    scs_free(a->c_aug);
    scs_free(a->tau);
    scs_free(a->qr_work);
    scs_free(a->jpvt);
    scs_free(a->W);
    scs_free(a->W_orig);
    scs_free(a->ipiv);
    scs_free(a->gamma_red);
    scs_free(a->c_top_save);
    scs_free(a->ir_res);
    scs_free(a->nrm_s_col);
    scs_free(a->nrm_y_col);
    scs_free(a->work);
    if (a->x_work) {
      scs_free(a->x_work);
    }
    scs_free(a);
  }
}

void aa_reset(AaWork *a) {
  /* Restore the logical state of a freshly calloc'd workspace.
   *
   * Most internal buffers are fully overwritten before they are read:
   *   - x, f, g_prev are re-seeded by init_accel_params on the next
   *     aa_apply call (which runs when iter == 0).
   *   - g, S/Y/D columns, A_aug/B_aug/c_aug, tau, jpvt, qr_work, W/W_orig,
   *     gamma_red, c_top_save, ir_res, work, ipiv, x_work are all
   *     rewritten inside update_accel_params / solve / relax each
   *     iteration before any read.
   *
   * The only buffers that require zeroing are nrm_{s,y}_col: they are
   * reduced over all `mem` slots in compute_regularization, and stale
   * entries from an earlier run would contaminate the Frobenius-norm
   * scale until every slot has been rewritten. */
  if (!a) {
    return;
  }
  if (a->verbosity > 0) {
    scs_printf("AA reset.\n");
  }
  a->iter = 0;
  a->success = 0;
  a->norm_g = 0;
  if (a->nrm_s_col) {
    memset(a->nrm_s_col, 0, sizeof(aa_float) * a->mem);
  }
  if (a->nrm_y_col) {
    memset(a->nrm_y_col, 0, sizeof(aa_float) * a->mem);
  }
}

AaStats aa_get_stats(const AaWork *a) {
  AaStats s;
  s.iter = a->iter;
  s.n_accept = a->n_accept;
  s.n_reject_lapack = a->n_reject_lapack;
  s.n_reject_rank0 = a->n_reject_rank0;
  s.n_reject_nonfinite = a->n_reject_nonfinite;
  s.n_reject_weight_cap = a->n_reject_weight_cap;
  s.n_safeguard_reject = a->n_safeguard_reject;
  s.last_rank = a->last_rank;
  s.last_aa_norm = a->last_aa_norm;
  s.last_regularization = a->last_regularization;
  return s;
}

#endif /* USE_LAPACK */
