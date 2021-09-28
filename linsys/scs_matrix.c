/* contains routines common to direct and indirect sparse solvers */
#include "scs_matrix.h"
#include "linalg.h"
#include "linsys.h"
#include "util.h"

#define MIN_NORMALIZATION_FACTOR (1e-4)
#define MAX_NORMALIZATION_FACTOR (1e4)
#define NUM_RUIZ_PASSES (25) /* additional passes don't help much */
#define NUM_L2_PASSES (1)    /* do one or zero, not more since not stable */

scs_int SCS(copy_matrix)(ScsMatrix **dstp, const ScsMatrix *src) {
  scs_int Anz = src->p[src->n];
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  if (!A) {
    return 0;
  }
  A->n = src->n;
  A->m = src->m;
  /* A values, size: NNZ A */
  A->x = (scs_float *)scs_malloc(sizeof(scs_float) * Anz);
  /* A row index, size: NNZ A */
  A->i = (scs_int *)scs_malloc(sizeof(scs_int) * Anz);
  /* A column pointer, size: n+1 */
  A->p = (scs_int *)scs_malloc(sizeof(scs_int) * (src->n + 1));
  if (!A->x || !A->i || !A->p) {
    return 0;
  }
  memcpy(A->x, src->x, sizeof(scs_float) * Anz);
  memcpy(A->i, src->i, sizeof(scs_int) * Anz);
  memcpy(A->p, src->p, sizeof(scs_int) * (src->n + 1));
  *dstp = A;
  return 1;
}

scs_int SCS(validate_lin_sys)(const ScsMatrix *A, const ScsMatrix *P) {
  scs_int i, j, r_max, Anz;
  if (!A->x || !A->i || !A->p) {
    scs_printf("data incompletely specified\n");
    return -1;
  }
  /* detects some errors in A col ptrs: */
  Anz = A->p[A->n];
  if (Anz > 0) {
    for (i = 0; i < A->n; ++i) {
      if (A->p[i] == A->p[i + 1]) {
        scs_printf("WARN: A->p (column pointers) not strictly increasing, "
                   "column %li empty\n",
                   (long)i);
      } else if (A->p[i] > A->p[i + 1]) {
        scs_printf("ERROR: A->p (column pointers) decreasing\n");
        return -1;
      }
    }
  }
  if (((scs_float)Anz / A->m > A->n) || (Anz < 0)) {
    scs_printf("Anz (nonzeros in A) = %li, outside of valid range\n",
               (long)Anz);
    return -1;
  }
  r_max = 0;
  for (i = 0; i < Anz; ++i) {
    if (A->i[i] > r_max) {
      r_max = A->i[i];
    }
  }
  if (r_max > A->m - 1) {
    scs_printf("number of rows in A inconsistent with input dimension\n");
    return -1;
  }
  if (P) {
    if (P->n != A->n) {
      scs_printf("P dimension = %li, inconsistent with n = %li\n", (long)P->n,
                 (long)A->n);
      return -1;
    }
    if (P->m != P->n) {
      scs_printf("P is not square\n");
      return -1;
    }
    for (j = 0; j < P->n; j++) { /* cols */
      for (i = P->p[j]; i < P->p[j + 1]; i++) {
        if (P->i[i] > j) { /* if row > */
          scs_printf("P is not upper triangular\n");
          return -1;
        }
      }
    }
  }
  return 0;
}

void SCS(free_scs_matrix)(ScsMatrix *A) {
  if (A) {
    scs_free(A->x);
    scs_free(A->i);
    scs_free(A->p);
    scs_free(A);
  }
}

static inline scs_float apply_limit(scs_float x) {
  /* need to bound to 1 for cols/rows of all zeros, otherwise blows up */
  x = x < MIN_NORMALIZATION_FACTOR ? 1.0 : x;
  x = x > MAX_NORMALIZATION_FACTOR ? MAX_NORMALIZATION_FACTOR : x;
  return x;
}

static void compute_ruiz_mats(ScsMatrix *P, ScsMatrix *A, scs_float *b,
                              scs_float *c, scs_float *Dt, scs_float *Et,
                              scs_float *s, scs_int *boundaries,
                              scs_int cone_boundaries_len) {
  scs_int i, j, kk, count, delta;
  scs_float wrk;

  /****************************  D  ****************************/

  /* initialize D */
  for (i = 0; i < A->m; ++i) {
    /* Dt[i] = 0.; */
    Dt[i] = ABS(b[i]);
  }

  /* calculate row norms */
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      Dt[A->i[j]] = MAX(Dt[A->i[j]], ABS(A->x[j]));
    }
  }

  /* accumulate D across each cone  */
  count = boundaries[0];
  for (i = 1; i < cone_boundaries_len; ++i) {
    delta = boundaries[i];
    wrk = SCS(norm_inf)(&(Dt[count]), delta);
    for (j = count; j < count + delta; ++j) {
      Dt[j] = wrk;
    }
    count += delta;
  }

  for (i = 0; i < A->m; ++i) {
    Dt[i] = SAFEDIV_POS(1.0, SQRTF(apply_limit(Dt[i])));
  }

  /****************************  E  ****************************/

  /* initialize E */
  for (i = 0; i < A->n; ++i) {
    /* Et[i] = 0.; */
    Et[i] = ABS(c[i]);
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

  /* calculate col norms, E */
  for (i = 0; i < A->n; ++i) {
    Et[i] = MAX(Et[i], SCS(norm_inf)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
    Et[i] = SAFEDIV_POS(1.0, SQRTF(apply_limit(Et[i])));
  }

  /* calculate s value */
  *s = MAX(SCS(norm_inf)(c, A->n), SCS(norm_inf)(b, A->m));
  *s = SAFEDIV_POS(1.0, SQRTF(apply_limit(*s)));
}

static void compute_l2_mats(ScsMatrix *P, ScsMatrix *A, scs_float *b,
                            scs_float *c, scs_float *Dt, scs_float *Et,
                            scs_float *s, scs_int *boundaries,
                            scs_int cone_boundaries_len) {
  scs_int i, j, kk, count, delta;
  scs_float wrk, norm_c, norm_b;

  /****************************  D  ****************************/

  /* initialize D */
  for (i = 0; i < A->m; ++i) {
    /* Dt[i] = 0.; */
    Dt[i] = b[i] * b[i];
  }

  /* calculate row norms */
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      Dt[A->i[j]] += A->x[j] * A->x[j];
    }
  }
  for (i = 0; i < A->m; ++i) {
    Dt[i] = SQRTF(Dt[i]); /* l2 norm of rows */
  }

  /* accumulate D across each cone  */
  count = boundaries[0];
  for (i = 1; i < cone_boundaries_len; ++i) {
    delta = boundaries[i];
    wrk = 0.;
    for (j = count; j < count + delta; ++j) {
      wrk += Dt[j];
    }
    wrk /= delta;
    for (j = count; j < count + delta; ++j) {
      Dt[j] = wrk;
    }
    count += delta;
  }

  for (i = 0; i < A->m; ++i) {
    Dt[i] = SAFEDIV_POS(1.0, SQRTF(apply_limit(Dt[i])));
  }

  /****************************  E  ****************************/

  /* initialize E */
  for (i = 0; i < A->n; ++i) {
    /* Et[i] = 0.; */
    Et[i] = c[i] * c[i];
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
        if (i != j) {
          Et[i] += wrk;
        }
      }
    }
  }

  /* calculate col norms, E */
  for (i = 0; i < A->n; ++i) {
    Et[i] += SCS(norm_sq)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]);
    Et[i] = SAFEDIV_POS(1.0, SQRTF(apply_limit(SQRTF(Et[i]))));
  }

  /* calculate s value */
  norm_c = SCS(norm_2)(c, A->n);
  norm_b = SCS(norm_2)(b, A->m);
  *s = SQRTF(norm_c * norm_c + norm_b * norm_b);
  *s = SAFEDIV_POS(1.0, SQRTF(apply_limit(*s)));
}

static void rescale(ScsMatrix *P, ScsMatrix *A, scs_float *b, scs_float *c,
                    scs_float *Dt, scs_float *Et, scs_float s, ScsScaling *scal,
                    scs_int *boundaries, scs_int cone_boundaries_len) {
  scs_int i, j;
  /* scale the rows of A with D */
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      A->x[j] *= Dt[A->i[j]];
    }
  }

  /* scale the cols of A with E */
  for (i = 0; i < A->n; ++i) {
    SCS(scale_array)(&(A->x[A->p[i]]), Et[i], A->p[i + 1] - A->p[i]);
  }

  if (P) {
    /* scale the rows of P with E */
    for (i = 0; i < P->n; ++i) {
      for (j = P->p[i]; j < P->p[i + 1]; ++j) {
        P->x[j] *= Et[P->i[j]];
      }
    }
    /* scale the cols of P with E */
    for (i = 0; i < P->n; ++i) {
      SCS(scale_array)(&(P->x[P->p[i]]), Et[i], P->p[i + 1] - P->p[i]);
    }
  }

  /* scale c */
  for (i = 0; i < A->n; ++i) {
    c[i] *= Et[i];
  }
  /* scale b */
  for (i = 0; i < A->m; ++i) {
    b[i] *= Dt[i];
  }

  /* Accumulate scaling */
  for (i = 0; i < A->m; ++i) {
    scal->D[i] *= Dt[i];
  }
  for (i = 0; i < A->n; ++i) {
    scal->E[i] *= Et[i];
  }

  /* Apply scaling */
  SCS(scale_array)(c, s, A->n);
  SCS(scale_array)(b, s, A->m);
  /* no need to scale P since primal_scale = dual_scale */
  /*
  if (P) {
    SCS(scale_array)(P->x, primal_scale, P->p[P->n]);
    SCS(scale_array)(P->x, 1.0 / dual_scale, P->p[P->n]);
  }
  */

  /* Accumulate scaling */
  scal->primal_scale *= s;
  scal->dual_scale *= s;
}

/* Will rescale as P -> EPE, A -> DAE, c -> sEc, b -> sDb, in-place.
 * Essentially trying to rescale this matrix:
 *
 * [P  A' c]   with   [E  0  0] on both sides (D, E diagonal)
 * [A  0  b]          [0  D  0]
 * [c' b' 0]          [0  0  s]
 *
 * which results in:
 *
 * [ EPE   EA'D  sEc ]
 * [ DAE    0    sDb ]
 * [ sc'E  sb'D   0  ]
 *
 * In other words D rescales the rows of A, b
 *                E rescales the cols of A and rows/cols of P, c'
 *
 * will repeatedly set: D^-1 ~ norm of rows of [ A  b ]
 *
 *                      E^-1 ~ norm of cols of [ P ]
 *                                             [ A ]
 *                                             [ c']
 *
 * `s` is incorporated into dual_scale and primal_scale
 *
 * The main complication is that D has to respect cone boundaries.
 *
 */
void SCS(normalize)(ScsMatrix *P, ScsMatrix *A, scs_float *b, scs_float *c,
                    ScsScaling *scal, scs_int *cone_boundaries,
                    scs_int cone_boundaries_len) {
  scs_int i;
  scs_float s;
  scs_float *Dt = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
  scs_float *Et = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
  scal->D = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
  scal->E = (scs_float *)scs_malloc(A->n * sizeof(scs_float));

#if VERBOSITY > 5
  SCS(timer) normalize_timer;
  SCS(tic)(&normalize_timer);
  scs_printf("normalizing A and P\n");
#endif

  /* init D, E */
  for (i = 0; i < A->m; ++i) {
    scal->D[i] = 1.;
  }
  for (i = 0; i < A->n; ++i) {
    scal->E[i] = 1.;
  }
  scal->primal_scale = 1.;
  scal->dual_scale = 1.;
  for (i = 0; i < NUM_RUIZ_PASSES; ++i) {
    compute_ruiz_mats(P, A, b, c, Dt, Et, &s, cone_boundaries,
                      cone_boundaries_len);
    rescale(P, A, b, c, Dt, Et, s, scal, cone_boundaries, cone_boundaries_len);
  }
  for (i = 0; i < NUM_L2_PASSES; ++i) {
    compute_l2_mats(P, A, b, c, Dt, Et, &s, cone_boundaries,
                    cone_boundaries_len);
    rescale(P, A, b, c, Dt, Et, s, scal, cone_boundaries, cone_boundaries_len);
  }
  scs_free(Dt);
  scs_free(Et);

#if VERBOSITY > 5
  scs_printf("finished normalizing A and P, time: %1.2es\n",
             SCS(tocq)(&normalize_timer) / 1e3);
  scs_printf("inf norm A %1.2e\n", SCS(norm_inf)(A->x, A->p[A->n]));
  if (P) {
    scs_printf("inf norm P %1.2e\n", SCS(norm_inf)(P->x, P->p[P->n]));
  }
  scs_printf("primal_scale %g\n", scal->primal_scale);
  scs_printf("dual_scale %g\n", scal->dual_scale);
  scs_printf("norm_b %g\n", SCS(norm_inf)(b, A->m));
  scs_printf("norm_c %g\n", SCS(norm_inf)(c, A->n));
  scs_printf("norm D %g\n", SCS(norm_inf)(scal->D, A->m));
  scs_printf("norm E %g\n", SCS(norm_inf)(scal->E, A->n));
#endif
}

void SCS(un_normalize)(ScsMatrix *A, ScsMatrix *P, const ScsScaling *scal) {
  scs_int i, j;
  scs_float *D = scal->D;
  scs_float *E = scal->E;
  for (i = 0; i < A->n; ++i) {
    SCS(scale_array)
    (&(A->x[A->p[i]]), 1. / E[i], A->p[i + 1] - A->p[i]);
  }
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      A->x[j] /= D[A->i[j]];
    }
  }
  if (P) {
    for (i = 0; i < P->n; ++i) {
      SCS(scale_array)
      (&(P->x[P->p[i]]), 1. / E[i], P->p[i + 1] - P->p[i]);
    }
    for (i = 0; i < P->n; ++i) {
      for (j = P->p[i]; j < P->p[i + 1]; ++j) {
        P->x[j] /= E[P->i[j]];
      }
    }
  }
}

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
  /* returns y += P x */
  scs_int p, j, i;
  scs_int n = P->n;
  scs_int *Pp = P->p;
  scs_int *Pi = P->i;
  scs_float *Px = P->x;
  /* y += P_upper x but skip diagonal entries*/
  for (j = 0; j < n; j++) { /* col */
    for (p = Pp[j]; p < Pp[j + 1]; p++) {
      i = Pi[p];    /* row */
      if (i != j) { /* skip the diagonal */
        y[i] += Px[p] * x[j];
      }
    }
  }
  /* y += P_lower x */
  SCS(accum_by_atrans)(P, x, y);
}
