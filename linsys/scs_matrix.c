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
  if (!src) {
    *dstp = SCS_NULL;
    return 1;
  }
  scs_int Anz = src->p[src->n];
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
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
  /* Disable this check which is slowish and typically just produces noise. */
  /*
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
  */
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

static void compute_ruiz_mats(ScsMatrix *P, ScsMatrix *A, scs_float *Dt,
                              scs_float *Et, ScsConeWork *cone) {
  scs_int i, j, kk;
  scs_float wrk;

  /****************************  D  ****************************/

  /* initialize D */
  for (i = 0; i < A->m; ++i) {
    Dt[i] = 0.;
    /* Dt[i] = ABS(b[i]); */
  }

  /* calculate row norms */
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      Dt[A->i[j]] = MAX(Dt[A->i[j]], ABS(A->x[j]));
    }
  }

  /* accumulate D across each cone  */
  SCS(enforce_cone_boundaries)(cone, Dt, &SCS(norm_inf));

  /* invert temporary vec to form D */
  for (i = 0; i < A->m; ++i) {
    Dt[i] = SAFEDIV_POS(1.0, SQRTF(apply_limit(Dt[i])));
  }

  /****************************  E  ****************************/

  /* initialize E */
  for (i = 0; i < A->n; ++i) {
    Et[i] = 0.;
    /* Et[i] = ABS(c[i]); */
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
}

static void compute_l2_mats(ScsMatrix *P, ScsMatrix *A, scs_float *Dt,
                            scs_float *Et, ScsConeWork *cone) {
  scs_int i, j, kk;
  scs_float wrk;

  /****************************  D  ****************************/

  /* initialize D */
  for (i = 0; i < A->m; ++i) {
    Dt[i] = 0.;
    /* Dt[i] = b[i] * b[i]; */
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
  SCS(enforce_cone_boundaries)(cone, Dt, &SCS(mean));

  for (i = 0; i < A->m; ++i) {
    Dt[i] = SAFEDIV_POS(1.0, SQRTF(apply_limit(Dt[i])));
  }

  /****************************  E  ****************************/

  /* initialize E */
  for (i = 0; i < A->n; ++i) {
    Et[i] = 0.;
    /* Et[i] = c[i] * c[i]; */
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
}

static void rescale(ScsMatrix *P, ScsMatrix *A, scs_float *Dt, scs_float *Et,
                    ScsScaling *scal, ScsConeWork *cone) {
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

  /* Accumulate scaling */
  for (i = 0; i < A->m; ++i) {
    scal->D[i] *= Dt[i];
  }
  for (i = 0; i < A->n; ++i) {
    scal->E[i] *= Et[i];
  }

  /* no need to scale P since later primal_scale = dual_scale */
  /*
  if (P) {
    SCS(scale_array)(P->x, primal_scale, P->p[P->n]);
    SCS(scale_array)(P->x, 1.0 / dual_scale, P->p[P->n]);
  }
  */
}

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
ScsScaling *SCS(normalize_a_p)(ScsMatrix *P, ScsMatrix *A, ScsConeWork *cone) {
  scs_int i;
  ScsScaling *scal = (ScsScaling *)scs_calloc(1, sizeof(ScsScaling));
  scs_float *Dt = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  scs_float *Et = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  scal->D = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  scal->E = (scs_float *)scs_calloc(A->n, sizeof(scs_float));

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
  for (i = 0; i < NUM_RUIZ_PASSES; ++i) {
    compute_ruiz_mats(P, A, Dt, Et, cone);
    rescale(P, A, Dt, Et, scal, cone);
  }
  for (i = 0; i < NUM_L2_PASSES; ++i) {
    compute_l2_mats(P, A, Dt, Et, cone);
    rescale(P, A, Dt, Et, scal, cone);
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
  scs_printf("norm D %g\n", SCS(norm_inf)(scal->D, A->m));
  scs_printf("norm E %g\n", SCS(norm_inf)(scal->E, A->n));
#endif
  return scal;
}

/*
void SCS(un_normalize_a_p)(ScsMatrix *A, ScsMatrix *P, const ScsScaling *scal) {
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
*/

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
