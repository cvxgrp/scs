#include "cones.h"
#include "linalg.h"
#include "scs_blas.h"
#include "scs_types.h"
#include "util.h" // just for timer

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements code for projecting onto the sum-of-largest eigenvalues
 * cone.
 *
 * Last modified: 25 August 2024.
 */

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(syev)(const char *jobz, const char *uplo, blas_int *n, scs_float *a,
                blas_int *lda, scs_float *w, scs_float *work, blas_int *lwork,
                blas_int *info);

void BLAS(scal)(const blas_int *n, const scs_float *sa, scs_float *sx,
                const blas_int *incx);

void BLAS(gemm)(const char *transa, const char *transb, blas_int *m,
                blas_int *n, blas_int *k, scs_float *alpha, scs_float *a,
                blas_int *lda, scs_float *b, blas_int *ldb, scs_float *beta,
                scs_float *c, blas_int *ldc);

#ifdef __cplusplus
}
#endif

// forward declaration
scs_int proj_sum_largest_cone_sorted(scs_float *t, scs_float *x, scs_int n,
                                     scs_int k);

void flip(scs_float *x, int n) {
  scs_float temp;
  for (int i = 0; i < n / 2; i++) {
    temp = x[i];
    x[i] = x[n - i - 1];
    x[n - i - 1] = temp;
  }
}

scs_int SCS(proj_sum_largest_evals)(scs_float *tX, scs_int n, scs_int k,
                                    ScsConeWork *c) {
  // tvX = [t, X], where X represents the lower triangular part of a matrix
  // stored in a compact form and off-diagonal elements have been scaled by
  // sqrt(2)
  scs_float *X = tX + 1;

  // ----------------------------------------------------------------------
  //                      compute eigendecomposition
  // ----------------------------------------------------------------------
  scs_int i;
  blas_int nb = (blas_int)n;
  blas_int nb_plus_one = (blas_int)(n + 1);
  blas_int one_int = 1;
  scs_float sqrt2 = sqrt(2.0);
  scs_float sqrt2_inv = 1.0 / sqrt2;
  scs_float *Xs = c->Xs;
  scs_float *e = c->e;
  scs_float *Z = c->Z;
  scs_float *work = c->work;
  blas_int lwork = c->lwork;
  blas_int info = 0;

  // copy lower triangular matrix into full matrix
  for (i = 0; i < n; ++i) {
    memcpy(&(Xs[i * (n + 1)]), &(X[i * n - ((i - 1) * i) / 2]),
           (n - i) * sizeof(scs_float));
  }

  // rescale diags by sqrt(2)
  BLAS(scal)(&nb, &sqrt2, Xs, &nb_plus_one);

  // Eigendecomposition. On exit, the lower triangular part of Xs stores
  // the eigenvectors. The vector e stores the eigenvalues in ascending
  // order (smallest eigenvalue first) */
  BLAS(syev)("Vectors", "Lower", &nb, Xs, &nb, e, work, &lwork, &info);
  if (info != 0) {
    scs_printf("WARN: LAPACK syev error, info = %i\n", (int)info);
    if (info < 0) {
      return info;
    }
  }

  // ----------------------------------------------------------------------
  //  Project onto spectral *vector* cone. Note that e is sqrt(2) times
  //  the eigenvalue vector we want to project. We therefore multiply
  //  tvX[0] by sqrt(2).
  // ----------------------------------------------------------------------
  tX[0] *= sqrt2;
  SPECTRAL_TIMING(SCS(timer) _timer; SCS(tic)(&_timer);)
  flip(e, n);
  scs_int status = proj_sum_largest_cone_sorted(&tX[0], e, n, k);
  flip(e, n);
  SPECTRAL_TIMING(c->tot_time_vec_cone_proj += SCS(tocq)(&_timer);)

  if (status < 0) {
    return status;
  }

  // ----------------------------------------------------------------------
  //             recover projection onto spectral *matrix* cone
  // ----------------------------------------------------------------------
  memcpy(c->work_sum_of_largest, Xs, n * n * sizeof(*Xs));
  for (i = 0; i < n; ++i) {
    BLAS(scal)(&nb, &e[i], &Xs[i * n], &one_int);
  }

  char transN = 'N';
  char transY = 'T';
  scs_float one = 1.0;
  scs_float zero = 0.0;

  // it is not safe to have overlapping matrices for dgemm_.
  BLAS(gemm)(&transN, &transY, &nb, &nb, &nb, &one, Xs, &nb,
             c->work_sum_of_largest, &nb, &zero, Z, &nb);

  BLAS(scal)(&nb, &sqrt2_inv, Z, &nb_plus_one);

  for (i = 0; i < n; ++i) {
    memcpy(&(X[i * n - ((i - 1) * i) / 2]), &(Z[i * (n + 1)]),
           (n - i) * sizeof(scs_float));
  }

  tX[0] *= sqrt2_inv;

  return 0;
}

