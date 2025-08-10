#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h"
#include "util.h"
#include "util_spectral_cones.h"

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 *
 * This file implements the projection onto the log-determinant cone.
 *
 * Last modified: 25 August 2024.
 */

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(syev)(const char *jobz, const char *uplo, blas_int *n, scs_float *a,
                blas_int *lda, scs_float *w, scs_float *work, blas_int *lwork,
                blas_int *info);
blas_int BLAS(syrk)(const char *uplo, const char *trans, const blas_int *n,
                    const blas_int *k, const scs_float *alpha,
                    const scs_float *a, const blas_int *lda,
                    const scs_float *beta, scs_float *c, const blas_int *ldc);
void BLAS(scal)(const blas_int *n, const scs_float *sa, scs_float *sx,
                const blas_int *incx);

#ifdef __cplusplus
}
#endif

// forward declare from log_cone_wrapper.c
scs_int log_cone_proj_wrapper(scs_float t0, scs_float v0, scs_float *x0,
                              scs_float *proj, scs_int n, scs_float *workspace,
                              Newton_stats *stats, bool *warm_start);

scs_int SCS(proj_logdet_cone)(scs_float *tvX, scs_int n, ScsConeWork *c,
                              scs_int offset, bool *warmstart) {
  Newton_stats *stats = &(c->newton_stats);

  // tvX = [t, v, X], where X represents the lower triangular part of a matrix
  // stored in a compact form and off-diagonal elements have been scaled by
  // sqrt(2)
  scs_float *X = tvX + 2;

#ifndef USE_LAPACK
  scs_printf("FAILURE: solving SDP but no blas/lapack libraries were found!\n");
  scs_printf("SCS will return nonsense!\n");
  SCS(scale_array)(X, NAN, n);
  return -1;
#endif

  // ----------------------------------------------------------------------
  //                      compute eigendecomposition
  // ----------------------------------------------------------------------
  scs_int i;
  blas_int nb = (blas_int)n;
  blas_int nb_plus_one = (blas_int)(n + 1);
  blas_int one_int = 1;
  scs_float zero = 0., one = 1.;
  scs_float sqrt2 = SQRTF(2.0);
  scs_float sqrt2_inv = 1.0 / sqrt2;
  scs_float *Xs = c->Xs;
  scs_float *Z = c->Z;
  scs_float *e = c->e;
  scs_float *work = c->work;
  blas_int lwork = c->lwork;
  blas_int info = 0;
  scs_float sq_eig_pos;

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
      scs_printf("entering LAPACK stuff!\n");
      return info;
    }
  }

  // ----------------------------------------------------------------------
  //  Project onto spectral *vector* cone. Note that e is sqrt(2) times
  //  the eigenvalue vector we want to project. We therefore multiply
  //  tvX[0] and tvX[1] by sqrt(2). The projection of sqrt(2) * (t0, v0,
  //  evals) is stored in current_log_proj
  // (or equivalently, in c->saved_log_projs + offset)
  // ----------------------------------------------------------------------
  scs_float *current_log_proj = c->saved_log_projs + offset;
  SPECTRAL_TIMING(SCS(timer) _timer; SCS(tic)(&_timer);)
  scs_int status =
      log_cone_proj_wrapper(sqrt2 * tvX[0], sqrt2 * tvX[1], e, current_log_proj,
                            n, c->work_logdet, stats, warmstart);
  SPECTRAL_TIMING(c->tot_time_vec_cone_proj += SCS(tocq)(&_timer);)

  if (status < 0) {
    return status;
  }
  // return immediately if the origin is the solution
  else if (status == IN_NEGATIVE_DUAL_CONE) {
    memset(tvX, 0, sizeof(scs_float) * ((n * (n + 1)) / 2 + 2));
    return 0;
  }

  // ----------------------------------------------------------------------
  //             recover projection onto spectral *matrix* cone
  // ----------------------------------------------------------------------
  scs_float *evals_proj = current_log_proj + 2;
  for (i = 0; i < n; ++i) {
    assert(evals_proj[i] >= 0);
    sq_eig_pos = SQRTF(evals_proj[i]);
    BLAS(scal)(&nb, &sq_eig_pos, &Xs[i * n], &one_int);
  }

  BLAS(syrk)("Lower", "NoTrans", &nb, &nb, &one, Xs, &nb, &zero, Z, &nb);
  BLAS(scal)(&nb, &sqrt2_inv, Z, &nb_plus_one);

  for (i = 0; i < n; ++i) {
    memcpy(&(X[i * n - ((i - 1) * i) / 2]), &(Z[i * (n + 1)]),
           (n - i) * sizeof(scs_float));
  }
  tvX[0] = sqrt2_inv * current_log_proj[0];
  tvX[1] = sqrt2_inv * current_log_proj[1];

  return 0;
}
