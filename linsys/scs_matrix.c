/* contains routines common to direct and indirect sparse solvers */

/* ======================== Includes / Constants ======================== */

#include "scs_matrix.h"
#include "cones.h"
#include "linalg.h"
#include "linsys.h"
#include "util.h"

#include <string.h>

#define MIN_NORMALIZATION_FACTOR (1e-4)
#define MAX_NORMALIZATION_FACTOR (1e4)
#define NUM_RUIZ_PASSES (25) /* additional passes don't help much */
#define NUM_L2_PASSES (1)    /* do one or zero, not more since not stable */

/* ======================== Matrix Copy / Free ======================== */

scs_int SCS(copy_matrix)(ScsMatrix **dstp, const ScsMatrix *src) {
  scs_int Anz;
  ScsMatrix *A;
  if (!src) {
    *dstp = SCS_NULL;
    return 1;
  }
  Anz = src->p[src->n];
  A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  if (!A) {
    return 0;
  }
  A->n = src->n;
  A->m = src->m;
  /* A values, size: NNZ A */
  A->x = (scs_float *)scs_calloc(Anz, sizeof(scs_float));
  /* A row index, size: NNZ A */
  A->i = (scs_int *)scs_calloc(Anz, sizeof(scs_int));
  /* A column pointer, size: n+1 */
  A->p = (scs_int *)scs_calloc(src->n + 1, sizeof(scs_int));
  if (!A->x || !A->i || !A->p) {
    scs_free(A->x);
    scs_free(A->i);
    scs_free(A->p);
    scs_free(A);
    return 0;
  }
  memcpy(A->x, src->x, sizeof(scs_float) * Anz);
  memcpy(A->i, src->i, sizeof(scs_int) * Anz);
  memcpy(A->p, src->p, sizeof(scs_int) * (src->n + 1));
  *dstp = A;
  return 1;
}

void SCS(free_scs_matrix)(ScsMatrix *A) {
  if (A) {
    scs_free(A->x);
    scs_free(A->i);
    scs_free(A->p);
    scs_free(A);
  }
}

/* ======================== Validation ======================== */

scs_int SCS(validate_lin_sys)(const ScsMatrix *A, const ScsMatrix *P) {
  scs_int i, j, Anz, Pnz;
  if (!A) {
    scs_printf("A matrix missing\n");
    return -1;
  }
  if (A->m <= 0 || A->n <= 0) {
    scs_printf("A matrix dimensions must be positive\n");
    return -1;
  }
  if (!A->x || !A->i || !A->p) {
    scs_printf("data incompletely specified\n");
    return -1;
  }
  if (A->p[0] != 0) {
    scs_printf("A->p[0] must equal 0\n");
    return -1;
  }
  for (j = 0; j < A->n; ++j) {
    if (A->p[j] < 0 || A->p[j] > A->p[j + 1]) {
      scs_printf("A->p (column pointers) must be nonnegative and "
                 "nondecreasing\n");
      return -1;
    }
  }
  Anz = A->p[A->n];
  if (((scs_float)Anz / A->m > A->n) || (Anz < 0)) {
    scs_printf("Anz (nonzeros in A) = %li, outside of valid range\n",
               (long)Anz);
    return -1;
  }
  for (i = 0; i < Anz; ++i) {
    if (A->i[i] < 0 || A->i[i] >= A->m) {
      scs_printf("A row index %li outside valid range [0, %li]\n",
                 (long)A->i[i], (long)A->m - 1);
      return -1;
    }
    if (!isfinite(A->x[i])) {
      scs_printf("A contains a non-finite entry\n");
      return -1;
    }
  }
  if (P) {
    if (!P->x || !P->i || !P->p) {
      scs_printf("P matrix incompletely specified\n");
      return -1;
    }
    if (P->n != A->n) {
      scs_printf("P dimension = %li, inconsistent with n = %li\n", (long)P->n,
                 (long)A->n);
      return -1;
    }
    if (P->m != P->n) {
      scs_printf("P is not square\n");
      return -1;
    }
    if (P->p[0] != 0) {
      scs_printf("P->p[0] must equal 0\n");
      return -1;
    }
    for (j = 0; j < P->n; ++j) {
      if (P->p[j] < 0 || P->p[j] > P->p[j + 1]) {
        scs_printf("P->p (column pointers) must be nonnegative and "
                   "nondecreasing\n");
        return -1;
      }
    }
    Pnz = P->p[P->n];
    if (((scs_float)Pnz / P->m > P->n) || (Pnz < 0)) {
      scs_printf("Pnz (nonzeros in P) = %li, outside of valid range\n",
                 (long)Pnz);
      return -1;
    }
    for (j = 0; j < P->n; j++) { /* cols */
      for (i = P->p[j]; i < P->p[j + 1]; i++) {
        if (P->i[i] < 0 || P->i[i] >= P->n) {
          scs_printf("P row index %li outside valid range [0, %li]\n",
                     (long)P->i[i], (long)P->n - 1);
          return -1;
        }
        if (P->i[i] > j) { /* if row > */
          scs_printf("P is not upper triangular\n");
          return -1;
        }
        if (!isfinite(P->x[i])) {
          scs_printf("P contains a non-finite entry\n");
          return -1;
        }
      }
    }
  }
  return 0;
}

/* ======================== Matrix-Vector Products ======================== */

void SCS(accum_by_atrans)(const ScsMatrix *A, const scs_float *x,
                          scs_float *y) {
  /* y += A'*x
     A in column compressed format
     parallelizes over columns (rows of A')
   */
  scs_int p, j;
  scs_int c1, c2;
  scs_float yj;
  scs_int n = A->n;
  scs_int *Ap = A->p;
  scs_int *Ai = A->i;
  scs_float *Ax = A->x;
#ifdef _OPENMP
#pragma omp parallel for private(p, c1, c2, yj)
#endif
  for (j = 0; j < n; j++) {
    yj = y[j];
    c1 = Ap[j];
    c2 = Ap[j + 1];
    for (p = c1; p < c2; p++) {
      yj += Ax[p] * x[Ai[p]];
    }
    y[j] = yj;
  }
}

void SCS(accum_by_a)(const ScsMatrix *A, const scs_float *x, scs_float *y) {
  /*y += A*x
    A in column compressed format
    */
  scs_int p, j, i;
  scs_int n = A->n;
  scs_int *Ap = A->p;
  scs_int *Ai = A->i;
  scs_float *Ax = A->x;
  for (j = 0; j < n; j++) { /* col */
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p]; /* row */
      y[i] += Ax[p] * x[j];
    }
  }
}

/* Since P is upper triangular need to be clever here */
void SCS(accum_by_p)(const ScsMatrix *P, const scs_float *x, scs_float *y) {
  /* returns y += P x where P is stored upper triangular (CSC).
   * Single pass: each stored entry (i,j) contributes to both y[i] (upper)
   * and y[j] (symmetric lower), halving NNZ traversals vs two-pass approach. */
  scs_int p, j;
  scs_int n = P->n;
  scs_int *Pp = P->p;
  scs_int *Pi = P->i;
  scs_float *Px = P->x;
  for (j = 0; j < n; j++) {
    for (p = Pp[j]; p < Pp[j + 1]; p++) {
      scs_int i = Pi[p];
      scs_float val = Px[p];
      y[i] += val * x[j]; /* upper triangle + diagonal */
      if (i != j) {
        y[j] += val * x[i]; /* symmetric lower triangle */
      }
    }
  }
}

/* ======================== Shared scs_solve_lin_sys ======================== */

scs_int SCS(lin_sys_shared_init)(ScsLinSysShared *sh, const ScsMatrix *A,
                                 const ScsMatrix *P, const scs_float *diag_r,
                                 scs_int cg, ScsLinSysPlainSolve solve,
                                 ScsLinSysResidual residual,
                                 ScsLinSysCorrect correct) {
  scs_int nm = MAX(A->n + A->m, 1);
  sh->A = A;
  sh->P = P;
  sh->diag_r = diag_r;
  sh->n = A->n;
  sh->m = A->m;
  sh->cg = cg;
  sh->solve = solve;
  sh->residual = residual;
  sh->correct = residual ? correct : SCS_NULL; /* reads the residual op's r */
  sh->b0 = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  sh->r = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  sh->c = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  if (!sh->b0 || !sh->r || !sh->c) {
    scs_printf("ERROR: linear solver workspace memory allocation failure\n");
    return -1;
  }
  return 0;
}

void SCS(lin_sys_shared_free)(ScsLinSysShared *sh) {
  if (sh) {
    scs_free(sh->b0);
    scs_free(sh->r);
    scs_free(sh->c);
    sh->b0 = sh->r = sh->c = SCS_NULL;
  }
}

/* r = b - K x with K = [R_x + P, A'; A, -R_y] from sh's A, P and diag_r
 * (SCS's products); returns ||r||_inf, and ||b||_inf in *nb unless nb is
 * SCS_NULL. r aliases neither b nor x. The statements and their order are
 * those of src/scs.c's per-iteration residual at 064f48a9 (iter_residual),
 * with the right-hand side read from b, so the values are the same to the
 * bit wherever that right-hand side was rounded before the subtraction (it
 * is a load here, so no compiler can fuse it into the subtraction). */
static scs_float kkt_residual(const ScsLinSysShared *sh, const scs_float *b,
                              const scs_float *x, scs_float *r,
                              scs_float *nb) {
  scs_int i, n = sh->n, m = sh->m;
  const scs_float *dr = sh->diag_r;
  scs_float rhs, nrhs = 0.;
  memset(r, 0, (n + m) * sizeof(scs_float));
  SCS(accum_by_atrans)(sh->A, &(x[n]), r);
  if (sh->P) {
    SCS(accum_by_p)(sh->P, x, r);
  }
  SCS(accum_by_a)(sh->A, x, &(r[n]));
  for (i = 0; i < n; ++i) {
    rhs = b[i];
    nrhs = MAX(nrhs, ABS(rhs));
    r[i] = rhs - (r[i] + dr[i] * x[i]);
  }
  for (i = n; i < n + m; ++i) {
    rhs = b[i];
    nrhs = MAX(nrhs, ABS(rhs));
    r[i] = rhs - (r[i] - dr[i] * x[i]);
  }
  if (nb) {
    *nb = nrhs;
  }
  return SCS(norm_inf)(r, n + m);
}

/* sh->r = b - K x (or where the solver keeps it, see ScsLinSysResidual) with
 * the solver's product if it has one, else SCS's; *nr = ||r||_inf, *nb =
 * ||b||_inf unless nb is SCS_NULL. Returns 0, or -1 if the solver's product
 * failed (with a message: that is an internal error, not a property of the
 * problem). */
static scs_int shared_residual(ScsLinSysShared *sh, ScsLinSysWork *w,
                               const scs_float *b, const scs_float *x,
                               scs_float *nr, scs_float *nb) {
  if (sh->residual) {
    if (sh->residual(w, b, x, sh->r, nr, nb) != 0) {
      scs_printf("ERROR: the linear solver's matrix product failed; please "
                 "open an issue at https://github.com/cvxgrp/scs/issues\n");
      return -1;
    }
    return 0;
  }
  *nr = kkt_residual(sh, b, x, sh->r, nb);
  return 0;
}

/* out = x + K^{-1} r for the x and r of the latest shared_residual call
 * (direct solvers): the solver's correct op, else the plain solve of a copy
 * of r plus x (one rounding: the scale is 1). out aliases neither x nor
 * sh->r. Returns the status of the solve. */
static scs_int shared_correct(ScsLinSysShared *sh, ScsLinSysWork *w,
                              const scs_float *x, scs_float *out,
                              scs_float tol) {
  scs_int status, nm = sh->n + sh->m;
  if (sh->correct) {
    return sh->correct(w, x, sh->r, out);
  }
  memcpy(out, sh->r, nm * sizeof(scs_float));
  status = sh->solve(w, out, SCS_NULL, tol);
  if (status == 0) {
    SCS(add_scaled_array)(out, x, nm, 1.);
  }
  return status;
}

/* scs_solve_lin_sys of linsys.h. b0 is filled only when a residual of the
 * solution will be needed (refinement or res). Every decision compares
 * residuals, so a call is deterministic whenever the solver's plain solve
 * and products are. */
scs_int SCS(lin_sys_shared_solve)(ScsLinSysShared *sh, ScsLinSysWork *w,
                                  scs_float *b, const scs_float *s,
                                  scs_float tol, scs_int max_refine_steps,
                                  scs_float *res) {
  scs_int nm = sh->n + sh->m, status, k, warm = 0, have_nr = 0, stop;
  scs_int need_b0 = max_refine_steps > 0 || res;
  scs_float nr = -1., nb = -1., nc;
  size_t vb = nm * sizeof(scs_float);
  if (sh->cg) {
    /* CG from s (its first n entries), or from zero, to tol: the solve of
     * SCS 3.3. No comparison with zero: a poor start costs CG iterations,
     * not accuracy. */
    if (need_b0) {
      memcpy(sh->b0, b, vb);
    }
    status = sh->solve(w, b, s, tol);
    if (status != 0) {
      return status;
    }
  } else {
    if (s) {
      /* correction form x = s + K^{-1}(b - K s), unless s is no better a
       * guess than zero (||b - K s|| >= ||b||): then b is untouched and the
       * plain solve below runs */
      status = shared_residual(sh, w, b, s, &nr, &nb);
      if (status != 0) {
        return status;
      }
      if (nr < nb) {
        if (need_b0) {
          memcpy(sh->b0, b, vb);
        }
        status = shared_correct(sh, w, s, b, tol);
        if (status != 0) {
          return status;
        }
        warm = 1;
      }
    }
    if (!warm) {
      if (need_b0) {
        memcpy(sh->b0, b, vb);
      }
      status = sh->solve(w, b, SCS_NULL, tol);
      if (status != 0) {
        return status;
      }
    }
  }
  if (max_refine_steps > 0) {
    status = shared_residual(sh, w, sh->b0, b, &nr, SCS_NULL);
    if (status != 0) {
      return status;
    }
    have_nr = 1;
    for (k = 0; k < max_refine_steps && nr > tol; ++k) {
      /* the candidate c: CG restarted on b from x, to tol (its residual is
       * recomputed below, since the one CG tracks can drift from it); or
       * x + K^{-1}(b - K x) with the residual of x just computed */
      if (sh->cg) {
        memcpy(sh->c, sh->b0, vb);
        status = sh->solve(w, sh->c, b, tol);
      } else {
        status = shared_correct(sh, w, b, sh->c, tol);
      }
      if (status != 0) {
        return status;
      }
      status = shared_residual(sh, w, sh->b0, sh->c, &nc, SCS_NULL);
      if (status != 0) {
        return status;
      }
      if (!(nc < nr)) {
        break; /* no better (or not finite): keep x */
      }
      memcpy(b, sh->c, vb);
      stop = nc > LIN_SYS_REFINE_CONTRACTION * nr; /* not halved */
      nr = nc;
      if (stop) {
        break;
      }
    }
  }
  if (res) {
    if (!have_nr) {
      status = shared_residual(sh, w, sh->b0, b, &nr, SCS_NULL);
      if (status != 0) {
        return status;
      }
    }
    *res = nr;
  }
  return 0;
}

/* ======================== Normalization Internals ======================== */

static inline scs_float apply_limit(scs_float x) {
  /* need to bound to 1 for cols/rows of all zeros, otherwise blows up */
  x = x < MIN_NORMALIZATION_FACTOR ? 1.0 : x;
  x = x > MAX_NORMALIZATION_FACTOR ? MAX_NORMALIZATION_FACTOR : x;
  return x;
}

/* Equilibrates the full homogeneous-embedding operator
 *   Q = [ 0   A'  c ]
 *       [-A   0   b ]
 *       [-c' -b'  0 ]
 * with the symmetric scaling diag(E, D, st): |Q| is symmetric so row and
 * column norms coincide, b enters the row norms, c the column norms, and
 * the scalar st on the tau slot replaces the old one-shot clamped sigma
 * heuristic in normalize_b_c. bt/ct are the running scaled copies of b/c
 * maintained by the caller; *st_out receives this pass's tau scaling. */
static void compute_ruiz_mats(ScsMatrix *P, ScsMatrix *A, const scs_float *bt,
                              const scs_float *ct, scs_float *Dt,
                              scs_float *Et, scs_float *st_out,
                              ScsConeWork *cone) {
  scs_int i, j, kk;
  scs_float wrk;

  /****************************  D  ****************************/

  /* initialize D with the tau-column contribution */
  for (i = 0; i < A->m; ++i) {
    Dt[i] = ABS(bt[i]);
  }

  /* calculate row norms */
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      Dt[A->i[j]] = MAX(Dt[A->i[j]], ABS(A->x[j]));
    }
  }

  /* accumulate D across each cone  */
  SCS(enforce_cone_boundaries)(cone, Dt, &SCS(norm_inf), 0);

  /* invert temporary vec to form D */
  for (i = 0; i < A->m; ++i) {
    Dt[i] = SQRTF(apply_limit(Dt[i]));
    Dt[i] = SAFEDIV_POS(1.0, Dt[i]);
  }

  /****************************  E  ****************************/

  /* initialize E with the tau-row contribution */
  for (i = 0; i < A->n; ++i) {
    Et[i] = ABS(ct[i]);
  }

  /* TODO: test not using P to determine scaling  */
  if (P) {
    /* compute norm of cols of P (symmetric upper triangular) */
    /* E = norm of cols of P */
    /* Compute maximum across columns */
    /* P(i, j) contributes to col j and col i (row i) due to symmetry */
    for (j = 0; j < P->n; j++) { /* cols */
      for (kk = P->p[j]; kk < P->p[j + 1]; kk++) {
        i = P->i[kk]; /* row */
        wrk = ABS(P->x[kk]);
        Et[j] = MAX(wrk, Et[j]);
        if (i != j) {
          Et[i] = MAX(wrk, Et[i]);
        }
      }
    }
  }

  /* calculate col norms, E — inline the norm_inf to avoid n BLAS call
   * overheads (Fortran ABI, pointer args, 1-based return) for short cols. */
  for (i = 0; i < A->n; ++i) {
    scs_float nm_a_col = 0.0, tmp;
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      tmp = ABS(A->x[j]);
      if (tmp > nm_a_col) nm_a_col = tmp;
    }
    Et[i] = MAX(Et[i], nm_a_col);
    Et[i] = SQRTF(apply_limit(Et[i]));
    Et[i] = SAFEDIV_POS(1.0, Et[i]);
  }

  /**************************  tau  ****************************/
  wrk = MAX(SCS(norm_inf)(bt, A->m), SCS(norm_inf)(ct, A->n));
  *st_out = SAFEDIV_POS(1.0, SQRTF(apply_limit(wrk)));
}

/* The l2 pass normalizes by RMS (root-mean-square over stored entries)
 * rather than the raw l2 norm. The raw l2 norm of a row/column grows like
 * sqrt(nnz) times its typical entry magnitude, so dividing by it crushes
 * dense or flat rows/columns by an extra nnz^(1/4) factor relative to
 * sparse ones -- a density artifact, not a magnitude signal (the Ruiz
 * linf passes are density-blind). RMS measures typical entry magnitude
 * while keeping the l2 pass's averaging robustness. Under the stacked
 * equilibration the b/c entries count toward their row/column averages,
 * and the tau norms average over their own lengths.
 */

static void compute_l2_mats(ScsMatrix *P, ScsMatrix *A, const scs_float *bt,
                            const scs_float *ct, scs_float *Dt, scs_float *Et,
                            scs_float *st_out, ScsConeWork *cone) {
  scs_int i, j, kk;
  scs_int rms = 1;
  scs_float wrk;
  scs_float *rcnt = SCS_NULL, *ecnt = SCS_NULL;
  rcnt = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  ecnt = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  if (!rcnt || !ecnt) {
    scs_free(rcnt);
    scs_free(ecnt);
    rcnt = ecnt = SCS_NULL;
    rms = 0; /* fall back to plain l2 on alloc failure */
  }

  /****************************  D  ****************************/

  /* initialize D with the tau-column contribution */
  for (i = 0; i < A->m; ++i) {
    Dt[i] = bt[i] * bt[i];
  }

  /* calculate row norms */
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      Dt[A->i[j]] += A->x[j] * A->x[j];
      if (rms) {
        rcnt[A->i[j]] += 1.;
      }
    }
  }
  for (i = 0; i < A->m; ++i) {
    if (rms) {
      Dt[i] /= (rcnt[i] + 1.); /* +1 for the bt (tau-column) entry */
    }
    Dt[i] = SQRTF(Dt[i]); /* l2 (or rms) norm of rows */
  }

  /* accumulate D across each cone  */
  SCS(enforce_cone_boundaries)(cone, Dt, &SCS(mean), 0);

  for (i = 0; i < A->m; ++i) {
    Dt[i] = SQRTF(apply_limit(Dt[i]));
    Dt[i] = SAFEDIV_POS(1.0, Dt[i]);
  }

  /****************************  E  ****************************/

  /* initialize E with the tau-row contribution */
  for (i = 0; i < A->n; ++i) {
    Et[i] = ct[i] * ct[i];
  }

  /* TODO: test not using P to determine scaling  */
  if (P) {
    /* compute norm of cols of P (symmetric upper triangular) */
    /* E = norm of cols of P */
    /* Compute maximum across columns */
    /* P(i, j) contributes to col j and col i (row i) due to symmetry */
    for (j = 0; j < P->n; j++) { /* cols */
      for (kk = P->p[j]; kk < P->p[j + 1]; kk++) {
        i = P->i[kk]; /* row */
        wrk = P->x[kk] * P->x[kk];
        Et[j] += wrk;
        if (rms) {
          ecnt[j] += 1.;
        }
        if (i != j) {
          Et[i] += wrk;
          if (rms) {
            ecnt[i] += 1.;
          }
        }
      }
    }
  }

  /* calculate col norms, E */
  for (i = 0; i < A->n; ++i) {
    Et[i] += SCS(norm_sq)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]);
    if (rms) {
      Et[i] /= (ecnt[i] + (scs_float)(A->p[i + 1] - A->p[i]) + 1.);
    }
    Et[i] = SQRTF(apply_limit(SQRTF(Et[i])));
    Et[i] = SAFEDIV_POS(1.0, Et[i]);
  }

  /**************************  tau  ****************************/
  if (rms) {
    wrk = MAX(SCS(norm_2)(bt, A->m) / SQRTF((scs_float)A->m),
              SCS(norm_2)(ct, A->n) / SQRTF((scs_float)A->n));
  } else {
    wrk = MAX(SCS(norm_2)(bt, A->m), SCS(norm_2)(ct, A->n));
  }
  *st_out = SAFEDIV_POS(1.0, SQRTF(apply_limit(wrk)));
  scs_free(rcnt);
  scs_free(ecnt);
}

static void rescale(ScsMatrix *P, ScsMatrix *A, scs_float *bt, scs_float *ct,
                    scs_float st, scs_float *Dt, scs_float *Et,
                    ScsScaling *scal, ScsConeWork *cone) {
  scs_int i, j;
  /* Fuse row and col scaling of A: A[i,j] *= Dt[i] * Et[j].
   * Single NNZ pass replaces two separate passes. */
  for (i = 0; i < A->n; ++i) {
    scs_float ei = Et[i];
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      A->x[j] *= Dt[A->i[j]] * ei;
    }
  }

  if (P) {
    /* Fuse row and col scaling of P: P[i,j] *= Et[i] * Et[j]. */
    for (i = 0; i < P->n; ++i) {
      scs_float ei = Et[i];
      for (j = P->p[i]; j < P->p[i + 1]; ++j) {
        P->x[j] *= Et[P->i[j]] * ei;
      }
    }
  }

  /* Scale the running b/c copies by their row/col factors and the tau
   * scalar (the b entries live at (row i, col tau), c at (row tau, col j)
   * of the stacked operator). */
  for (i = 0; i < A->m; ++i) {
    bt[i] *= Dt[i] * st;
  }
  for (i = 0; i < A->n; ++i) {
    ct[i] *= Et[i] * st;
  }

  /* Accumulate scaling */
  for (i = 0; i < A->m; ++i) {
    scal->D[i] *= Dt[i];
  }
  for (i = 0; i < A->n; ++i) {
    scal->E[i] *= Et[i];
  }
  scal->tau_scale *= st;

  /* no need to scale P since later primal_scale = dual_scale */
  /*
  if (P) {
    SCS(scale_array)(P->x, primal_scale, P->p[P->n]);
    SCS(scale_array)(P->x, 1.0 / dual_scale, P->p[P->n]);
  }
  */
}

/* ======================== Normalization Public API ======================== */

/* Will rescale as P -> EPE, A -> DAE in-place.
 * Essentially trying to rescale this matrix:
 *
 * [P  A']   with   [E  0 ] on both sides (D, E diagonal)
 * [A  0 ]          [0  D ]
 *
 * which results in:
 *
 * [ EPE   EA'D ]
 * [ DAE    0   ]
 *
 * In other words D rescales the rows of A
 *                E rescales the cols of A and rows/cols of P
 *
 * will repeatedly set: D^-1 ~ norm of rows of [ A ]
 *
 *                      E^-1 ~ norm of cols of [ P ]
 *                                             [ A ]
 *
 * The main complication is that D has to respect cone boundaries.
 *
 */
ScsScaling *SCS(normalize_a_p)(ScsMatrix *P, ScsMatrix *A, const scs_float *b,
                               const scs_float *c, ScsConeWork *cone) {
  scs_int i;
  scs_float st;
  ScsScaling *scal = (ScsScaling *)scs_calloc(1, sizeof(ScsScaling));
  scs_float *Dt = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  scs_float *Et = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  /* running scaled copies of b, c (originals must not be modified here;
   * their actual scaling happens per-solve in normalize_b_c) */
  scs_float *bt = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
  scs_float *ct = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
  if (!scal || !Dt || !Et || !bt || !ct) {
    scs_free(scal);
    scs_free(Dt);
    scs_free(Et);
    scs_free(bt);
    scs_free(ct);
    return SCS_NULL;
  }
  memcpy(bt, b, A->m * sizeof(scs_float));
  memcpy(ct, c, A->n * sizeof(scs_float));
  scal->D = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  scal->E = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  if (!scal->D || !scal->E) {
    scs_free(scal->D);
    scs_free(scal->E);
    scs_free(scal);
    scs_free(Dt);
    scs_free(Et);
    return SCS_NULL;
  }

#if VERBOSITY > 5
  SCS(timer) normalize_timer;
  SCS(tic)(&normalize_timer);
  scs_printf("normalizing A and P\n");
#endif

  /* init D, E */
  scal->m = A->m;
  for (i = 0; i < A->m; ++i) {
    scal->D[i] = 1.;
  }
  scal->n = A->n;
  for (i = 0; i < A->n; ++i) {
    scal->E[i] = 1.;
  }
  scal->primal_scale = 1.;
  scal->dual_scale = 1.;
  scal->tau_scale = 1.;
  for (i = 0; i < NUM_RUIZ_PASSES; ++i) {
    compute_ruiz_mats(P, A, bt, ct, Dt, Et, &st, cone);
    rescale(P, A, bt, ct, st, Dt, Et, scal, cone);
  }
  for (i = 0; i < NUM_L2_PASSES; ++i) {
    compute_l2_mats(P, A, bt, ct, Dt, Et, &st, cone);
    rescale(P, A, bt, ct, st, Dt, Et, scal, cone);
  }

  scs_free(Dt);
  scs_free(Et);
  scs_free(bt);
  scs_free(ct);

#if VERBOSITY > 5
  scs_printf("finished normalizing A and P, time: %1.2es\n",
             SCS(tocq)(&normalize_timer) / 1e3);
  scs_printf("inf norm A %1.2e\n", SCS(norm_inf)(A->x, A->p[A->n]));
  if (P) {
    scs_printf("inf norm P %1.2e\n", SCS(norm_inf)(P->x, P->p[P->n]));
  }
  scs_printf("primal_scale %g\n", scal->primal_scale);
  scs_printf("dual_scale %g\n", scal->dual_scale);
  scs_printf("norm D %g\n", SCS(norm_inf)(scal->D, A->m));
  scs_printf("norm E %g\n", SCS(norm_inf)(scal->E, A->n));
#endif
  return scal;
}
