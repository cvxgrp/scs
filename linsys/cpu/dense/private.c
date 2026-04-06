#include "private.h"

#include <string.h>

/* LAPACK prototypes */
#ifdef __cplusplus
extern "C" {
#endif

void BLAS(getrf)(blas_int *m, blas_int *n, scs_float *a, blas_int *lda,
                 blas_int *ipiv, blas_int *info);
void BLAS(getrs)(char *trans, blas_int *n, blas_int *nrhs, const scs_float *a,
                 blas_int *lda, const blas_int *ipiv, scs_float *b,
                 blas_int *ldb, blas_int *info);

#ifdef __cplusplus
}
#endif

const char *scs_get_lin_sys_method(void) {
  return "dense-direct-dgetrf";
}

/* Fill the dense KKT matrix from sparse A, P and diagonal R:
 *
 *   KKT = [ (R_x + P)   A' ]
 *         [    A        -R_y ]
 *
 * Stored column-major. Uses full (unsymmetric) matrix for dgetrf.
 */
static void form_dense_kkt(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i, j, k;
  scs_int n = p->n;
  scs_int m = p->m;
  scs_int nm = p->n_plus_m;
  scs_float *K = p->kkt;
  const ScsMatrix *A = p->A;
  const ScsMatrix *P = p->P;

  memset(K, 0, (size_t)nm * nm * sizeof(scs_float));

  /* Top-left block: R_x + P (n x n) */
  /* Diagonal: R_x */
  for (i = 0; i < n; ++i) {
    K[i * nm + i] = diag_r[i];
    p->diag_p[i] = 0.;
  }
  /* Add P (symmetric, stored upper triangular CSC) */
  if (P) {
    for (j = 0; j < n; ++j) {
      for (k = P->p[j]; k < P->p[j + 1]; ++k) {
        i = P->i[k];
        K[j * nm + i] += P->x[k]; /* upper: (i, j) */
        if (i == j) {
          p->diag_p[j] = P->x[k];
        } else {
          K[i * nm + j] += P->x[k]; /* lower: (j, i) symmetric fill */
        }
      }
    }
  }

  /* Off-diagonal blocks: A (m x n) in bottom-left, A' (n x m) in top-right */
  /* A is stored CSC: column j has entries A->i[k], A->x[k] */
  for (j = 0; j < n; ++j) {
    for (k = A->p[j]; k < A->p[j + 1]; ++k) {
      i = A->i[k];
      K[j * nm + (n + i)] = A->x[k]; /* bottom-left: A at row n+i, col j */
      K[(n + i) * nm + j] = A->x[k]; /* top-right: A' at row j, col n+i */
    }
  }

  /* Bottom-right block: -R_y (m x m) diagonal */
  for (i = 0; i < m; ++i) {
    K[(n + i) * nm + (n + i)] = -diag_r[n + i];
  }
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  blas_int info, bn;
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p)
    return SCS_NULL;

  p->n = A->n;
  p->m = A->m;
  p->n_plus_m = A->n + A->m;
  p->A = A;
  p->P = P;

  p->kkt = (scs_float *)scs_calloc((size_t)p->n_plus_m * p->n_plus_m,
                                    sizeof(scs_float));
  p->ipiv = (blas_int *)scs_calloc(p->n_plus_m, sizeof(blas_int));
  p->diag_p = (scs_float *)scs_calloc(p->n, sizeof(scs_float));

  if (!p->kkt || !p->ipiv || !p->diag_p) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Form and factorize */
  form_dense_kkt(p, diag_r);
  bn = (blas_int)p->n_plus_m;
  BLAS(getrf)(&bn, &bn, p->kkt, &bn, p->ipiv, &info);
  if (info != 0) {
    scs_printf("Error in dense LU factorization (dgetrf), info = %d\n",
               (int)info);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  return p;
}

/* Solve KKT * x = b, solution stored in b */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  blas_int bn = (blas_int)p->n_plus_m;
  blas_int one = 1;
  blas_int info;
  char trans = 'N';
  BLAS(getrs)(&trans, &bn, &one, p->kkt, &bn, p->ipiv, b, &bn, &info);
  if (info != 0) {
    scs_printf("Error in dense solve (dgetrs), info = %d\n", (int)info);
  }
  return (scs_int)info;
}

/* Update diagonal R entries and re-factorize */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  blas_int info, bn;
  form_dense_kkt(p, diag_r);
  bn = (blas_int)p->n_plus_m;
  BLAS(getrf)(&bn, &bn, p->kkt, &bn, p->ipiv, &info);
  if (info != 0) {
    scs_printf("Error in dense LU re-factorization (dgetrf), info = %d\n",
               (int)info);
    return (scs_int)info;
  }
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    scs_free(p->kkt);
    scs_free(p->ipiv);
    scs_free(p->diag_p);
    scs_free(p);
  }
}
