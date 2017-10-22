#include "scs.h"
#include "accel.h"
#include "linalg.h"
#include "scs_blas.h"
#include "util.h"

/* This file uses Anderson acceleration to improve the convergence of the
 * ADMM iteration z^+ = \phi(z). At each iteration we need to solve a (small)
 * linear system, we do this using LAPACK, first forming the normal equations
 * and using ?posv (fastest, but bad numerical stability), if that fails we
 * switch to using ?gels, which uses a QR factorization (slower, but better
 * numerically). If this fails then we just don't do any acceleration this
 * iteration, however we could fall back further to ?gelsy or other more
 * robust methods if we wanted to.
 */

#define REGULARIZATION (0.0)

struct SCS_ACCEL_WORK {
#ifdef LAPACK_LIB_FOUND
  scs_float *d_f;
  scs_float *d_g;
  scs_float *d_x;
  scs_float *f;
  scs_float *g;
  scs_float *x;
  scs_float *sol;
  scs_float *scratch;
  scs_float *mat;
  blas_int * ipiv;
  scs_int k, l;
#endif
  scs_float total_accel_time;
};

#ifdef LAPACK_LIB_FOUND
void BLAS(gemv)(const char *trans, const blas_int *m, const blas_int *n,
    const scs_float *alpha, const scs_float *a, const blas_int *lda,
    const scs_float *x, const blas_int *incx, const scs_float *beta,
    scs_float *y, const blas_int *incy);
void BLAS(gemm)(const char *transa, const char *transb, blas_int *m, blas_int *
    n, blas_int *k, scs_float *alpha, scs_float *a, blas_int *lda,
    scs_float *b, blas_int *ldb, scs_float *beta, scs_float *c, blas_int
    *ldc);
void BLAS(gesv) (blas_int * n, blas_int * nrhs, scs_float * a,
    blas_int * lda, blas_int * ipiv, scs_float * b, blas_int * ldb, blas_int * info);

void update_accel_params(ScsWork *w, scs_int idx) {
  DEBUG_FUNC
  scs_float *d_f = w->accel->d_f;
  scs_float *d_g = w->accel->d_g;
  scs_float *d_x = w->accel->d_x;
  scs_float *f = w->accel->f;
  scs_float *g = w->accel->g;
  scs_float *x = w->accel->x;
  scs_int l = w->m + w->n + 1;
  /* copy g_prev into idx col of d_g */
  memcpy(&(d_g[idx * 2 * l]), g, sizeof(scs_float) * 2 * l);
  /* copy f_prev into idx col of d_f */
  memcpy(&(d_f[idx * 2 * l]), f, sizeof(scs_float) * 2 * l);
  /* copy x_prev into idx col of d_x */
  memcpy(&(d_x[idx * 2 * l]), x, sizeof(scs_float) * 2 * l);
  /* x = [u_prev;v_prev] */
  memcpy(x, w->u_prev, sizeof(scs_float) * l);
  memcpy(&(x[l]), w->v_prev, sizeof(scs_float) * l);
  /* g = [u;v] */
  memcpy(g, w->u, sizeof(scs_float) * l);
  memcpy(&(g[l]), w->v, sizeof(scs_float) * l);
  /* calculcate f = g - x */
  memcpy(f, g, sizeof(scs_float) * 2 * l);
  add_scaled_array(f, x, 2 * l, -1.0);
  /* idx col of d_g = g_prev - g */
  add_scaled_array(&(d_g[idx * 2 * l]), g, 2 * l, -1);
  /* idx col of d_f = f_prev - f */
  add_scaled_array(&(d_f[idx * 2 * l]), f, 2 * l, -1);
  /* idx col of d_x = x_prev - x */
  add_scaled_array(&(d_x[idx * 2 * l]), x, 2 * l, -1);
  RETURN;
}

ScsAccelWork *init_accel(ScsWork *w) {
  DEBUG_FUNC
  ScsAccelWork *a = scs_calloc(1, sizeof(ScsAccelWork));
  if (!a) {
    RETURN SCS_NULL;
  }
  a->l = w->m + w->n + 1;
  /* k = lookback - 1 since we use the difference form
     of anderson acceleration, and so there is one fewer var in lin sys.
     Use MIN to prevent not full rank matrices */
  a->k = MIN(w->n, w->stgs->acceleration_lookback - 1);
  if (a->k <= 0) {
    RETURN a;
  }
  a->d_f = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->d_g = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->d_x = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->f = scs_calloc(2 * a->l, sizeof(scs_float));
  a->g = scs_calloc(2 * a->l, sizeof(scs_float));
  a->x = scs_calloc(2 * a->l, sizeof(scs_float));
  a->mat = scs_calloc(a->k * a->k, sizeof(scs_float));
  a->sol = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->scratch = scs_malloc(sizeof(scs_float) * a->k);
  a->ipiv = scs_malloc(sizeof(blas_int) * a->k);
  a->total_accel_time = 0.0;
  if (!a->d_f || !a->d_g || !a->f || !a->g || !a->scratch || !a->sol || !a->d_x || !a->x || !a->scratch || !a->ipiv || !a->mat) {
    free_accel(a);
    a = SCS_NULL;
  }
  RETURN a;
}

scs_int solve_with_gesv(ScsAccelWork * a, scs_int len) {
  DEBUG_FUNC
  scs_int i;
  blas_int info;
  blas_int twol = 2 * a->l;
  blas_int one = 1;
  blas_int blen = (blas_int) len;
  scs_float neg_onef = -1.0;
  scs_float onef = 1.0;
  scs_float zerof = 0.0;
  scs_float * d_x = a->d_x;
  scs_float * d_f = a->d_f;
  scs_float regularization = REGULARIZATION;
  if (regularization > 0.) {
    memset(a->mat, 0, len * len * sizeof(scs_float));
    for (i = 0; i < len; ++i) {
      a->mat[i * len + i] = 1.0;
    }
  }
  /* mat = dX'*dF */
  BLAS(gemm)("Trans", "NoTrans", &blen, &blen, &twol, &onef, d_x, &twol, d_f, &twol, &regularization, a->mat, &blen);
  /* scratch = dX' f */
  BLAS(gemv)("Trans", &twol, &blen, &onef, d_x, &twol, a->f, &one, &zerof, a->scratch, &one);
  /* scratch = (dX'dF) \ dX' f */
  BLAS(gesv)(&blen, &one, a->mat, &blen, a->ipiv, a->scratch, &blen, &info);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * a->l);
  /* sol = sol - dG * scratch */
  BLAS(gemv)("NoTrans", &twol, &blen, &neg_onef, a->d_g, &twol, a->scratch, &one, &onef, a->sol, &one);
  RETURN (scs_int) info;
}

scs_int accelerate(ScsWork *w, scs_int iter) {
  DEBUG_FUNC
  scs_int l = w->accel->l;
  scs_int k = w->accel->k;
  scs_int info;
  timer accel_timer;
  tic(&accel_timer);
  if (k <= 0) {
    RETURN 0;
  }
  /* update df, d_g, d_x, f, g, x */
  update_accel_params(w, iter % k);
  if (iter == 0) {
    RETURN 0;
  }
  /* solve linear system, new point stored in sol */
  info = solve_with_gesv(w->accel, MIN(iter, k));
  /* set [u;v] = sol */
  memcpy(w->u, w->accel->sol, sizeof(scs_float) * l);
  memcpy(w->v, &(w->accel->sol[l]), sizeof(scs_float) * l);
  w->accel->total_accel_time += tocq(&accel_timer);
  /* add check that info == 0 and fallback */
  RETURN info;
}

void free_accel(ScsAccelWork *a) {
  DEBUG_FUNC
  if (a) {
    if (a->d_f) scs_free(a->d_f);
    if (a->d_g) scs_free(a->d_g);
    if (a->d_x) scs_free(a->d_x);
    if (a->f) scs_free(a->f);
    if (a->g) scs_free(a->g);
    if (a->x) scs_free(a->x);
    if (a->sol) scs_free(a->sol);
    if (a->scratch) scs_free(a->scratch);
    if (a->mat) scs_free(a->mat);
    if (a->ipiv) scs_free(a->ipiv);
    scs_free(a);
  }
  RETURN;
}

#else

ScsAccelWork *init_accel(ScsWork *w) {
  ScsAccelWork *a = scs_malloc(sizeof(ScsAccelWork));
  a->total_accel_time = 0.0;
  RETURN a;
}

void free_accel(ScsAccelWork *a) {
  if (a) {
    scs_free(a);
  }
}

scs_int accelerate(ScsWork *w, scs_int iter) { RETURN 0; }
#endif

char *get_accel_summary(const ScsInfo *info, ScsAccelWork *a) {
  DEBUG_FUNC
  char *str = scs_malloc(sizeof(char) * 64);
  sprintf(str, "\tAcceleration: avg step time: %1.2es\n",
          a->total_accel_time / (info->iter + 1) / 1e3);
  a->total_accel_time = 0.0;
  RETURN str;
}
