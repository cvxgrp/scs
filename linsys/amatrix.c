/* contains routines common to direct and indirect sparse solvers */
#include "amatrix.h"

#include "linsys.h"

#define MIN_SCALE (1e-4)
#define MAX_SCALE (1e4)
#define NUM_SCALE_PASSES 10 /* additional passes don't help much */

/* Typically l2 equilibration works better than l_inf (Ruiz) */
/* Though more experimentation is needed */
/* #define RUIZ 1 */

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
        scs_printf(
            "WARN: A->p (column pointers) not strictly increasing, "
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

#if EXTRA_VERBOSE > 0
static void print_matrix(const ScsMatrix *A) {
  scs_int i, j;
  /* TODO: this limit of 250 is to prevent clogging stdout */
  if (A->p[A->n] < 250) {
    scs_printf("\n");
    for (i = 0; i < A->n; ++i) {
      scs_printf("Col %li: ", (long)i);
      for (j = A->p[i]; j < A->p[i + 1]; j++) {
        scs_printf("A[%li,%li] = %4f, ", (long)A->i[j], (long)i, A->x[j]);
      }
      scs_printf("norm col = %4f\n",
                 SCS(norm)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
    }
    scs_printf("norm A = %4f\n", SCS(norm)(A->x, A->p[A->n]));
  }
}
#endif

static inline scs_float apply_limit(scs_float x) {
  x = x < MIN_SCALE ? 1.0 : x;
  x = x > MAX_SCALE ? MAX_SCALE : x;
  return x;
}

/* Ruiz style rescaling using inf norms */
void SCS(_normalize)(ScsMatrix *A, ScsMatrix *P, const ScsSettings *stgs,
                     const ScsCone *k, ScsScaling *scal) {
  scs_float *D = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
  scs_float *E = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
  scs_float *Dt = (scs_float *)scs_malloc(A->m * sizeof(scs_float));
  scs_float *Et = (scs_float *)scs_malloc(A->n * sizeof(scs_float));
  scs_int i, j, kk, l, count, delta, *boundaries;
  scs_int num_boundaries = SCS(get_cone_boundaries)(k, &boundaries);
  scs_float wrk, norm_a, norm_p;

#if EXTRA_VERBOSE > 0
  SCS(timer) normalize_timer;
  SCS(tic)(&normalize_timer);
  scs_printf("normalizing A and P\n");
  print_matrix(A);
  if (P) print_matrix(P);
#endif

/* Will rescale as P -> EPE, A -> DAE.
 * Essentially trying to rescale this matrix:
 *
 * [P  A' c]   with   [E  0  0] on both sides (D, E diagonal)
 * [A  0  b]          [0  D  0]
 * [c' b' 0]          [0  0  1]
 *
 * which results in:
 *
 * [EPE EA'D  Ec]
 * [DAE  0    Db]
 * [c'E  b'D   0]
 *
 * (Though we don't use b and c to pick D, E so that we can reuse the same
 * factorization for many b,c combinations).
 *
 * In other words D rescales the rows of A
 *                E rescales the cols of A and rows/cols of P
 *
 * will repeatedly set D^-1 ~ norm of rows of A
 *                     E^-1 ~ norm of cols of [P]
 *                                            [A]
 *
 * The main complication is that D has to respect cone boundaries.
 *
 * */

/* Balance A and P to begin */
#ifdef RUIZ
  norm_a = SCS(norm_inf)(A->x, A->p[A->n]);
#else
  norm_a = SCS(norm)(A->x, A->p[A->n]); /* should be square to approx A'A ? */
#endif
  norm_a = apply_limit(norm_a);
  scal->scale_p = 0.;
  if (P) {
#ifdef RUIZ
    norm_p = SCS(norm_inf)(P->x, P->p[P->n]);
#else
    norm_p = SCS(norm)(P->x, P->p[P->n]);
#endif
    /* If P is set but is zero then same as not set */
    if (norm_p > 0.) {
      norm_p = apply_limit(norm_p);
      scal->scale_p = norm_a / norm_p;
      SCS(scale_array)(P->x, scal->scale_p, P->p[P->n]);
    }
  }

  for (l = 0; l < NUM_SCALE_PASSES; ++l) {
    memset(D, 0, A->m * sizeof(scs_float));
    memset(E, 0, A->n * sizeof(scs_float));
    /* calculate row norms */
    for (i = 0; i < A->n; ++i) {
      for (j = A->p[i]; j < A->p[i + 1]; ++j) {
#ifdef RUIZ
        D[A->i[j]] = MAX(D[A->i[j]], ABS(A->x[j]));
#else
        /* l2 normalization row norms sqrt(n / m) */
        D[A->i[j]] += A->m * A->x[j] * A->x[j] / A->n;
#endif
      }
    }

    /* accumulate D across each cone  */
    count = boundaries[0];
    for (i = 1; i < num_boundaries; ++i) {
      delta = boundaries[i];
#ifdef RUIZ
      wrk = SCS(norm_inf)(&(D[count]), delta);
#else
      wrk = 0;
      for (j = count; j < count + delta; ++j) {
        wrk += D[j];
      }
      wrk /= delta;
#endif
      for (j = count; j < count + delta; ++j) {
        D[j] = wrk;
      }
      count += delta;
    }

    for (i = 0; i < A->m; ++i) {
#ifdef RUIZ
      D[i] = apply_limit(SQRTF(D[i]));
#else
      D[i] = apply_limit(SQRTF(SQRTF(D[i])));
#endif
    }

    if (P) {
      /* compute inf norm of cols of P (symmetric upper triangular) */
      /* E = inf norm of cols of P */
      /* Compute maximum across columns */
      /* P(i, j) contributes to col j and col i (row i) due to symmetry */
      for (j = 0; j < P->n; j++) { /* cols */
        for (kk = P->p[j]; kk < P->p[j + 1]; kk++) {
          i = P->i[kk]; /* row */
          wrk = ABS(P->x[kk]);
#ifdef RUIZ
          E[j] = MAX(wrk, E[j]);
          if (i != j) {
            E[i] = MAX(wrk, E[i]);
          }
#else
          E[j] += wrk * wrk;
          if (i != j) {
            E[i] += wrk * wrk;
          }
#endif
        }
      }
    }

    /* calculate col norms, E */
    for (i = 0; i < A->n; ++i) {
#ifdef RUIZ
      E[i] = MAX(E[i], SCS(norm_inf)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
      E[i] = apply_limit(SQRTF(E[i]));
#else
      E[i] += SCS(norm_sq)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]);
      E[i] = apply_limit(SQRTF(SQRTF(E[i])));
#endif
    }

    /* scale the rows of A with D */
    for (i = 0; i < A->n; ++i) {
      for (j = A->p[i]; j < A->p[i + 1]; ++j) {
        A->x[j] /= D[A->i[j]];
      }
    }

    /* scale the cols of A with E */
    for (i = 0; i < A->n; ++i) {
      SCS(scale_array)(&(A->x[A->p[i]]), 1.0 / E[i], A->p[i + 1] - A->p[i]);
    }

    if (P) {
      /* scale the rows of P with E */
      for (i = 0; i < P->n; ++i) {
        for (j = P->p[i]; j < P->p[i + 1]; ++j) {
          P->x[j] /= E[P->i[j]];
        }
      }
      /* scale the cols of P with E */
      for (i = 0; i < P->n; ++i) {
        SCS(scale_array)(&(P->x[P->p[i]]), 1.0 / E[i], P->p[i + 1] - P->p[i]);
      }
    }

    /* Accumulate scaling */
    for (i = 0; i < A->m; ++i) {
      Dt[i] = (l == 0) ? D[i] : Dt[i] * D[i];
    }
    for (i = 0; i < A->n; ++i) {
      Et[i] = (l == 0) ? E[i] : Et[i] * E[i];
    }
  }
  scs_free(boundaries);
  scs_free(D);
  scs_free(E);

#ifdef RUIZ
  scal->norm_a = SCS(norm_inf)(A->x, A->p[A->n]);
#else
  scal->norm_a = SCS(norm)(A->x, A->p[A->n]);
#endif
  /* scale up by d->SCALE if not equal to 1 */
  if (stgs->scale != 1) {
    SCS(scale_array)(A->x, stgs->scale, A->p[A->n]);
    if (P) {
      SCS(scale_array)(P->x, stgs->scale, P->p[P->n]);
    }
  }

  scal->D = Dt;
  scal->E = Et;

#if EXTRA_VERBOSE > 0
  scs_printf("finished normalizing A and P, time: %1.2es\n",
             SCS(tocq)(&normalize_timer) / 1e3);
  print_matrix(A);
  scs_printf("inf norm A %1.2e\n", SCS(norm_inf)(A->x, A->p[A->n]));
  if (P) {
    print_matrix(P);
    scs_printf("inf norm P %1.2e\n", SCS(norm_inf)(P->x, P->p[P->n]));
  }
#endif
}

void SCS(_un_normalize)(ScsMatrix *A, ScsMatrix *P, const ScsSettings *stgs,
                        const ScsScaling *scal) {
  scs_int i, j;
  scs_float *D = scal->D;
  scs_float *E = scal->E;
  for (i = 0; i < A->n; ++i) {
    SCS(scale_array)
    (&(A->x[A->p[i]]), E[i] / stgs->scale, A->p[i + 1] - A->p[i]);
  }
  for (i = 0; i < A->n; ++i) {
    for (j = A->p[i]; j < A->p[i + 1]; ++j) {
      A->x[j] *= D[A->i[j]];
    }
  }
  if (P) {
    for (i = 0; i < P->n; ++i) {
      SCS(scale_array)
      (&(P->x[P->p[i]]), E[i] / stgs->scale, P->p[i + 1] - P->p[i]);
    }
    for (i = 0; i < P->n; ++i) {
      for (j = P->p[i]; j < P->p[i + 1]; ++j) {
        P->x[j] *= E[P->i[j]];
      }
    }
  }
}

scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n) {
  scs_int i, nz = 0;
  scs_float nz2 = 0;
  if (!p || !c) {
    return (-1);
  } /* check inputs */
  for (i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i]; /* also in scs_float to avoid scs_int overflow */
    c[i] = p[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p[n] = nz;
  return nz2; /* return sum (c [0..n-1]) */
}

void SCS(_accum_by_atrans)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                           const scs_float *x, scs_float *y) {
  /* y += A'*x
     A in column compressed format
     parallelizes over columns (rows of A')
   */
  scs_int p, j;
  scs_int c1, c2;
  scs_float yj;
#if EXTRA_VERBOSE > 0
  SCS(timer) mult_by_atrans_timer;
  SCS(tic)(&mult_by_atrans_timer);
#endif
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
#if EXTRA_VERBOSE > 0
  scs_printf("trans mat mul time: %1.2es\n",
             SCS(tocq)(&mult_by_atrans_timer) / 1e3);
#endif
}

void SCS(_accum_by_a)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                      const scs_float *x, scs_float *y, scs_int skip_diag) {
  /*y += A*x
    A in column compressed format
   */
  scs_int p, j, i;
#if EXTRA_VERBOSE > 0
  SCS(timer) mult_by_a_timer;
  SCS(tic)(&mult_by_a_timer);
#endif
  /* have skip_diag out here for compiler optimization */
  if (skip_diag) {
    for (j = 0; j < n; j++) { /* col */
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        i = Ai[p]; /* row */
        if (i != j) {
          y[i] += Ax[p] * x[j];
        }
      }
    }
  } else {
    for (j = 0; j < n; j++) { /* col */
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        i = Ai[p]; /* row */
        y[i] += Ax[p] * x[j];
      }
    }
  }
#if EXTRA_VERBOSE > 0
  scs_printf("mat mult time: %1.2es\n", SCS(tocq)(&mult_by_a_timer) / 1e3);
#endif
}

/* Since P is upper triangular need to be clever here */
void SCS(accum_by_p)(const ScsMatrix *P, ScsLinSysWork *p, const scs_float *x,
                     scs_float *y) {
  /* returns y += P x */
  /* y += P_upper x but skip diagonal entries*/
  SCS(_accum_by_a)(P->n, P->x, P->i, P->p, x, y, 1);
  /* y += P_lower x */
  SCS(_accum_by_atrans)(P->n, P->x, P->i, P->p, x, y);
}

