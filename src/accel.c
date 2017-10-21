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
  scs_float *Q;
  scs_float *R;
  scs_float *delta;
  scs_float *dummy_row;
  scs_float *mat;
  scs_int k, l;
#endif
  scs_float total_accel_time;
};

#ifdef LAPACK_LIB_FOUND
void BLAS(gemv)(const char *trans, const blas_int *m, const blas_int *n,
    const scs_float *alpha, const scs_float *a, const blas_int *lda,
    const scs_float *x, const blas_int *incx, const scs_float *beta,
    scs_float *y, const blas_int *incy);
void BLAS(geqrf)(blas_int *m, blas_int *n, scs_float *a, blas_int *lda,
    scs_float *tau, scs_float *work, blas_int *lwork,
    blas_int *info);
void BLAS(orgqr)(blas_int *m, blas_int *n, blas_int *k, scs_float *a, blas_int *lda,
    scs_float *tau, scs_float *work, blas_int *lwork,
    blas_int *info);
void BLAS(trsv)(const char *uplo, const char *trans, const char *diag,
    blas_int *n, scs_float *a, blas_int *lda, scs_float *x,
    blas_int *incx);
void BLAS(rotg)(scs_float *a, scs_float *b, scs_float *c, scs_float *s);
void BLAS(rot)(const blas_int *n, scs_float *x, const blas_int *incx,
    scs_float *y, const blas_int *incy, const scs_float *c,
    const scs_float *s);
scs_float BLAS(nrm2)(const blas_int *n, scs_float *x, const blas_int *incx);
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
  scs_float *delta = w->accel->delta;
  scs_int l = w->m + w->n + 1;
  /* copy old col into delta */
  memcpy(delta, &(d_f[idx * 2 * l]), sizeof(scs_float) * 2 * l);
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
  /* delta = delta - new col */
  add_scaled_array(delta, &(d_f[idx * 2 * l]), 2 * l, -1.0);
  /* delta = new - old */
  scale_array(delta, -1.0, 2 * l);
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
  a->Q = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->R = scs_calloc(a->k * a->k, sizeof(scs_float));
  a->mat = scs_calloc(a->k * a->k, sizeof(scs_float));
  a->dummy_row = scs_calloc(a->k, sizeof(scs_float));
  a->sol = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->scratch = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->delta = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->total_accel_time = 0.0;
  if (!a->d_f || !a->d_g || !a->f || !a->g || !a->scratch || !a->sol ||
      !a->delta || !a->R) {
    free_accel(a);
    a = SCS_NULL;
  }
  RETURN a;
}

void qrfactorize(ScsAccelWork *a) {
  DEBUG_FUNC
  scs_int l = a->l;
  scs_int i;
  blas_int twol = 2 * l;
  blas_int bk = (blas_int)a->k;
  blas_int neg_one = -1;
  blas_int info;
  blas_int lwork;
  scs_float worksize;
  scs_float *work;
  scs_float *tau = scs_malloc(a->k * sizeof(scs_float));
  scs_float *Q = a->Q;
  memcpy(Q, a->d_f, sizeof(scs_float) * a->k * 2 * a->l);
  BLAS(geqrf)(&twol, &bk, Q, &twol, tau, &worksize, &neg_one, &info);
  lwork = (blas_int)worksize;
  work = scs_malloc(lwork * sizeof(scs_float));
  BLAS(geqrf)(&twol, &bk, Q, &twol, tau, work, &lwork, &info);
  scs_free(work);
  for (i = 0; i < a->k; ++i) {
    memcpy(&(a->R[i * a->k]), &(Q[i * a->l * 2]), sizeof(scs_float) * (i + 1));
  }
  BLAS(orgqr)(&twol, &bk, &bk, Q, &twol, tau, &worksize, &neg_one, &info);
  lwork = (blas_int)worksize;
  work = scs_malloc(lwork * sizeof(scs_float));
  BLAS(orgqr)(&twol, &bk, &bk, Q, &twol, tau, work, &lwork, &info);
  scs_free(work);
  scs_free(tau);
  RETURN;
}

void update_factorization(ScsAccelWork *a, scs_int idx) {
  DEBUG_FUNC
  scs_float *Q = a->Q;
  scs_float *R = a->R;
  scs_float *u = a->delta;
  scs_float *w = a->scratch;
  scs_float *dummy_row = a->dummy_row;
  blas_int one = 1;
  blas_int bk = (blas_int)a->k;
  blas_int twol = (blas_int)2 * a->l;
  scs_float zerof = 0.0;
  scs_float onef = 1.0;
  scs_float neg_onef = -1.0;
  scs_float nrm_u;
  scs_float c;
  scs_float s;

  scs_int k = a->k;
  scs_int l = a->l;
  scs_float r1, r2;
  scs_int i, ridx;

  memset(dummy_row, 0, k * sizeof(scs_float));
  /* w = Q' * delta, size k: col of R */
  BLAS(gemv)("Trans", &twol, &bk, &onef, a->Q, &twol, u, &one, &zerof, w, &one);
  /* u = delta - Q * w = d_f * R^-1 w, size m: col of Q */
  BLAS(gemv)("no_trans", &twol, &bk, &neg_onef, a->Q, &twol, w, &one, &onef, u, &one);
  /* nrm_u = ||u|| */
  nrm_u = BLAS(nrm2)(&twol, u, &one);
  /* u = u / ||u|| */
  scale_array(u, 1.0 / nrm_u, 2 * l);
  /* R col += w */
  add_scaled_array(&(R[idx * k]), w, k, 1.0);

  /* Givens rotations, start with fake bottom row of R, extra col of Q */
  ridx = k * idx + k - 1;
  r1 = R[ridx];
  r2 = nrm_u;
  dummy_row[idx] = nrm_u;
  BLAS(rotg)(&r1, &r2, &c, &s);
  BLAS(rot)(&bk, &(R[k - 1]), &bk, dummy_row, &one, &c, &s);
  BLAS(rot)(&twol, &(Q[2 * l * (k - 1)]), &one, u, &one, &c, &s);

  /* Walk up the spike, R finishes upper Hessenberg */
  for (i = k; i > idx + 1; --i) {
    ridx = k * idx + i - 1;
    r1 = R[ridx - 1];
    r2 = R[ridx];
    BLAS(rotg)(&(r1), &(r2), &c, &s);
    BLAS(rot)(&bk, &(R[i - 2]), &bk, &(R[i - 1]), &bk, &c, &s);
    BLAS(rot)
    (&twol, &(Q[2 * l * (i - 2)]), &one, &(Q[2 * l * (i - 1)]), &one, &c, &s);
  }

  /* Walk down the sub-diagonal, R finishes upper triangular */
  for (i = idx + 1; i < k - 1; ++i) {
    ridx = k * i + i;
    r1 = R[ridx];
    r2 = R[ridx + 1];
    BLAS(rotg)(&r1, &r2, &c, &s);
    BLAS(rot)(&bk, &(R[i]), &bk, &(R[i + 1]), &bk, &c, &s);
    BLAS(rot)
    (&twol, &(Q[2 * l * i]), &one, &(Q[2 * l * (i + 1)]), &one, &c, &s);
  }

  /* Finish fake bottom row of R, extra col of Q */
  BLAS(rotg)(&(R[k * k - 1]), &(dummy_row[k - 1]), &c, &s);
  BLAS(rot)(&twol, &(Q[2 * l * (k - 1)]), &one, u, &one, &c, &s);
  RETURN;
}

scs_int solve_with_gesv(ScsAccelWork * a, scs_int len) {
  DEBUG_FUNC
  scs_int i;
  blas_int info;
  blas_int twol = 2 * a->l;
  blas_int one = 1;
  blas_int blen = (blas_int) len;
  scs_float negOnef = -1.0;
  scs_float onef = 1.0;
  scs_float zerof = 0.0;
  scs_float * d_x = a->d_x;
  scs_float * d_f = a->d_f;
  //memset(a->mat, 0, len * len * sizeof(scs_float));
  /* 
  if (reg > 0.) {
    for (i = 0; i < len; ++i) {
      a->mat[i * len + i] = 1.;
    }
  }
  */
  /* mat = dX'*dF */
  /* a->d_x = a->d_f; */
  BLAS(gemm)("Trans", "NoTrans", &blen, &blen, &twol, &onef, d_x, &twol, d_f, &twol, &zerof, a->mat, &blen);
  /* scratch = dX' f */
  BLAS(gemv)("Trans", &twol, &blen, &onef, d_x, &twol, a->f, &one, &zerof, a->scratch, &one);
  /* scratch = (dX'dF) \ dX' f */
  blas_int * ipiv = scs_malloc(sizeof(blas_int) * blen);
  BLAS(gesv)(&blen, &one, a->mat, &blen, ipiv, a->scratch, &blen, &info);
  scs_free(ipiv);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * a->l);
  /* sol = sol - dG * scratch */
  BLAS(gemv)("NoTrans", &twol, &blen, &negOnef, a->d_g, &twol, a->scratch, &one, &onef, a->sol, &one);
  RETURN (scs_int) info;
}

scs_int accelerate(ScsWork *w, scs_int iter) {
  DEBUG_FUNC
  scs_int l = w->accel->l;
  scs_int k = w->accel->k;
  scs_float *sol;
  scs_int idx, info;
  timer accel_timer;
  if (k <= 0) {
    RETURN 0;
  }
  tic(&accel_timer);
  idx = k - 1 - iter % k;
  update_accel_params(w, idx);
  /* iter < k doesn't do any acceleration until iters hit k */
  if (iter < k + 1) {
    /* if (iter == 0) { */
    RETURN 0;
  }
  /* solve linear system, new point stored in sol */
  info = solve_with_gesv(w->accel, k);
  /* set [u;v] = sol */
  memcpy(w->u, w->accel->sol, sizeof(scs_float) * l);
  memcpy(w->v, &(w->accel->sol[l]), sizeof(scs_float) * l);
  w->accel->total_accel_time += tocq(&accel_timer);
  RETURN 0;
}

void free_accel(ScsAccelWork *a) {
  DEBUG_FUNC
  if (a) {
    if (a->d_f) scs_free(a->d_f);
    if (a->d_g) scs_free(a->d_g);
    if (a->d_x) scs_free(a->d_g);
    if (a->f) scs_free(a->f);
    if (a->g) scs_free(a->g);
    if (a->x) scs_free(a->x);
    if (a->sol) scs_free(a->sol);
    if (a->scratch) scs_free(a->scratch);
    if (a->Q) scs_free(a->Q);
    if (a->R) scs_free(a->R);
    if (a->delta) scs_free(a->delta);
    if (a->dummy_row) scs_free(a->dummy_row);
    if (a->mat) scs_free(a->mat);
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
