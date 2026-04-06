/*
 * Dense direct linear system solver using the Gram matrix reduction:
 *
 *   [R_x + P   A'] [x]   [rx]
 *   [  A     -R_y] [y] = [ry]
 *
 * is reduced to:
 *
 *   G x = rx + A' R_y^{-1} ry,    where G = R_x + P + A' R_y^{-1} A
 *   y = R_y^{-1} (A x - ry)
 *
 * G is n x n symmetric positive definite, factorized via Cholesky (dpotrf).
 * A is stored as a dense m x n matrix for fast dgemv operations.
 */

#include "private.h"

#include <string.h>

/* LAPACK / BLAS prototypes */
#ifdef __cplusplus
extern "C" {
#endif

/* Cholesky factorization */
void BLAS(potrf)(const char *uplo, blas_int *n, scs_float *a, blas_int *lda,
                 blas_int *info);
/* Cholesky solve */
void BLAS(potrs)(const char *uplo, blas_int *n, blas_int *nrhs,
                 const scs_float *a, blas_int *lda, scs_float *b,
                 blas_int *ldb, blas_int *info);
/* Matrix-vector multiply: y = alpha * op(A) * x + beta * y */
void BLAS(gemv)(const char *trans, blas_int *m, blas_int *n, scs_float *alpha,
                const scs_float *a, blas_int *lda, const scs_float *x,
                blas_int *incx, scs_float *beta, scs_float *y, blas_int *incy);
/* Symmetric rank-k update: C = alpha * A' * A + beta * C */
void BLAS(syrk)(const char *uplo, const char *trans, blas_int *n, blas_int *k,
                scs_float *alpha, const scs_float *a, blas_int *lda,
                scs_float *beta, scs_float *c, blas_int *ldc);

#ifdef __cplusplus
}
#endif

const char *scs_get_lin_sys_method(void) {
  return "dense-direct-cholesky";
}

/* Convert sparse CSC matrix A (m x n) to dense column-major array */
static void sparse_to_dense(const ScsMatrix *A, scs_float *A_dense) {
  scs_int j, k;
  scs_int m = A->m;
  scs_int n = A->n;
  memset(A_dense, 0, (size_t)m * n * sizeof(scs_float));
  for (j = 0; j < n; ++j) {
    for (k = A->p[j]; k < A->p[j + 1]; ++k) {
      A_dense[j * m + A->i[k]] = A->x[k];
    }
  }
}

/* Form Gram matrix G = R_x + P + A' diag(r_y_inv) A (upper triangle).
 * Uses BLAS dsyrk for the A' diag(r_y_inv) A product. */
static void form_gram(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i, j, k;
  scs_int n = p->n;
  scs_int m = p->m;
  scs_float *G = p->G;
  const ScsMatrix *P = p->P;

  /* Compute r_y_inv = 1 / R_y and form scaled A: tmp_m used as scratch */
  for (i = 0; i < m; ++i) {
    p->r_y_inv[i] = 1.0 / diag_r[n + i];
  }

  /* Form S = diag(sqrt(r_y_inv)) * A in pre-allocated workspace.
   * S is m x n, stored column-major. */
  scs_float *S = p->S;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      S[j * m + i] = SQRTF(p->r_y_inv[i]) * p->A_dense[j * m + i];
    }
  }

  /* G = S' * S = A' diag(r_y_inv) A  (upper triangle) */
  {
    blas_int bn = (blas_int)n;
    blas_int bm = (blas_int)m;
    scs_float one = 1.0;
    scs_float zero = 0.0;
    char uplo = 'U';
    char trans = 'T';
    BLAS(syrk)(&uplo, &trans, &bn, &bm, &one, S, &bm, &zero, G, &bn);
  }

  /* Add R_x to diagonal */
  for (i = 0; i < n; ++i) {
    G[i * n + i] += diag_r[i];
    p->diag_p[i] = 0.;
  }

  /* Add P (symmetric, stored upper triangular CSC) */
  if (P) {
    for (j = 0; j < n; ++j) {
      for (k = P->p[j]; k < P->p[j + 1]; ++k) {
        i = P->i[k];
        G[j * n + i] += P->x[k]; /* upper triangle: i <= j */
        if (i == j) {
          p->diag_p[j] = P->x[k];
        }
      }
    }
  }
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  blas_int info, bn;
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p) {
    return SCS_NULL;
  }

  p->n = A->n;
  p->m = A->m;
  p->P = P;

  p->A_dense = (scs_float *)scs_calloc((size_t)A->m * A->n, sizeof(scs_float));
  p->G = (scs_float *)scs_calloc((size_t)A->n * A->n, sizeof(scs_float));
  p->r_y_inv = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  p->tmp_m = (scs_float *)scs_calloc(A->m, sizeof(scs_float));
  p->S = (scs_float *)scs_calloc((size_t)A->m * A->n, sizeof(scs_float));
  p->diag_p = (scs_float *)scs_calloc(A->n, sizeof(scs_float));

  if (!p->A_dense || !p->G || !p->r_y_inv || !p->tmp_m || !p->S || !p->diag_p) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Convert sparse A to dense once */
  sparse_to_dense(A, p->A_dense);

  /* Form and factorize Gram matrix */
  form_gram(p, diag_r);
  bn = (blas_int)p->n;
  {
    char uplo = 'U';
    BLAS(potrf)(&uplo, &bn, p->G, &bn, &info);
  }
  if (info != 0) {
    scs_printf("Error in dense Cholesky factorization (dpotrf), info = %d\n",
               (int)info);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  return p;
}

/*
 * Solve:
 *   x = G^{-1} (rx + A' R_y^{-1} ry)
 *   y = R_y^{-1} (A x - ry)
 *
 * Input b = [rx; ry], solution stored in b = [x; y].
 */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  blas_int bn = (blas_int)p->n;
  blas_int bm = (blas_int)p->m;
  blas_int one = 1;
  blas_int info;
  scs_float alpha, beta;
  scs_int i;

  scs_float *rx = b;
  scs_float *ry = b + p->n;

  /* tmp_m = R_y^{-1} ry */
  for (i = 0; i < p->m; ++i) {
    p->tmp_m[i] = p->r_y_inv[i] * ry[i];
  }

  /* rx += A' * tmp_m   (rx = rx + A' R_y^{-1} ry) */
  {
    char trans = 'T';
    alpha = 1.0;
    beta = 1.0;
    BLAS(gemv)(&trans, &bm, &bn, &alpha, p->A_dense, &bm, p->tmp_m, &one,
               &beta, rx, &one);
  }

  /* Solve G x = rx, result in rx (= b[:n]) */
  {
    char uplo = 'U';
    BLAS(potrs)(&uplo, &bn, &one, p->G, &bn, rx, &bn, &info);
  }
  if (info != 0) {
    scs_printf("Error in dense Cholesky solve (dpotrs), info = %d\n",
               (int)info);
    return (scs_int)info;
  }

  /* ry = A x - ry */
  {
    char trans = 'N';
    alpha = 1.0;
    beta = -1.0;
    BLAS(gemv)(&trans, &bm, &bn, &alpha, p->A_dense, &bm, rx, &one, &beta, ry,
               &one);
  }

  /* ry = R_y^{-1} (A x - ry) = y */
  for (i = 0; i < p->m; ++i) {
    ry[i] *= p->r_y_inv[i];
  }

  return 0;
}

/* Update diagonal R entries: re-form Gram matrix and re-factorize */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  blas_int info, bn;
  form_gram(p, diag_r);
  bn = (blas_int)p->n;
  {
    char uplo = 'U';
    BLAS(potrf)(&uplo, &bn, p->G, &bn, &info);
  }
  if (info != 0) {
    scs_printf(
        "Error in dense Cholesky re-factorization (dpotrf), info = %d\n",
        (int)info);
    return (scs_int)info;
  }
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    scs_free(p->A_dense);
    scs_free(p->G);
    scs_free(p->r_y_inv);
    scs_free(p->tmp_m);
    scs_free(p->S);
    scs_free(p->diag_p);
    scs_free(p);
  }
}
