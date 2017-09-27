#include "accel.h"
#include "scs.h"
#include "scs_blas.h"

/* Not clear if this should just be 0. */
#define ACCEL_REGULARIZATION (0.)

/* This file uses Anderson acceleration to improve the convergence of the
 * ADMM iteration z^+ = \phi(z). At each iteration we need to solve a (small)
 * linear system, we do this using LAPACK, first forming the normal equations
 * and using ?posv (fastest, but bad numerical stability), if that fails we
 * switch to using ?gels, which uses a QR factorization (slower, but better
 * numerically). If this fails then we just don't do any acceleration this
 * iteration, however we could fall back further to ?gelsy or other more
 * robust methods if we wanted to.
 */

struct SCS_ACCEL {
#ifdef LAPACK_LIB_FOUND
  scs_float *dF;
  scs_float *dG;
  scs_float *f;
  scs_float *g;
  scs_float *sol;
  scs_float *theta;
  scs_float *X; /* dF'dF for posv, or QR of dF for gels */
  scs_float *wrk; /* gels only: workspace for QR decomp */
  blasint worksize;
  scs_int k, l, bad_numerics;
#endif
  scs_float totalAccelTime;
};

#ifdef LAPACK_LIB_FOUND
void BLAS(gemv)(const char *trans, const blasint *m, const blasint *n,
                const scs_float *alpha, const scs_float *a, const blasint *lda,
                const scs_float *x, const blasint *incx, const scs_float *beta,
                scs_float *y, const blasint *incy);
void BLAS(syrk)(const char *uplo, const char *trans, blasint *n, blasint *k,
                scs_float *alpha, scs_float *a, blasint *lda, scs_float *beta,
                scs_float *c, blasint *ldc);
void BLAS(posv) (const char *uplo, blasint * n, blasint * nrhs, scs_float * a,
                 blasint * lda, scs_float * b, blasint * ldb, blasint * info);
void BLAS(gels)(const char *trans, const blasint *m, const blasint *n,
                const blasint *nrhs, scs_float *a, const blasint *lda,
                scs_float *b, const blasint *ldb, scs_float *work,
                const blasint *lwork, blasint *info);

scs_int solve_with_gels(Accel * a, scs_int len) {
  DEBUG_FUNC
  blasint info;
  scs_int l = a->l;
  blasint twol = 2 * l;
  blasint one = 1;
  blasint blen = (blasint) len;
  scs_float negOnef = -1.0;
  scs_float onef = 1.0;
  memcpy(a->theta, a->f, sizeof(scs_float) * 2 * l);
  memcpy(a->X, a->dF, sizeof(scs_float) * 2 * l * len);
  /* solve dF theta = f */
  BLAS(gels)("NoTrans", &twol, &blen, &one, a->X, &twol, a->theta, &twol, a->wrk, &(a->worksize), &info);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * l);
  /* sol = sol - dG * theta */
  BLAS(gemv)("NoTrans", &twol, &blen, &negOnef, a->dG, &twol, a->theta, &one, &onef, a->sol, &one);
  RETURN info;
}

scs_int solve_with_posv(Accel * a, scs_int len) {
  DEBUG_FUNC
  scs_int i;
  blasint info;
  blasint twol = 2 * a->l;
  blasint one = 1;
  blasint blen = (blasint) len;
  scs_float negOnef = -1.0;
  scs_float onef = 1.0;
  scs_float zerof = 0.0;
  scs_float reg = ACCEL_REGULARIZATION;
  memset(a->X, 0, len * len * sizeof(scs_float));
  if (reg > 0.) {
    for (i = 0; i < len; ++i) {
      a->X[i * len + i] = 1.;
    }
  }
  /* X = dF'*dF */
  BLAS(syrk)("Lower", "Transpose", &blen, &twol, &onef, a->dF, &twol, &reg, a->X, &blen);
  /* theta = dF' f */
  BLAS(gemv)("Trans", &twol, &blen, &onef, a->dF, &twol, a->f, &one, &zerof, a->theta, &one);
  /* theta = (dF'dF) \ dF' f */
  BLAS(posv)("Lower", &blen, &one, a->X, &blen, a->theta, &blen, &info);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * a->l);
  /* sol = sol - dG * theta */
  BLAS(gemv)("NoTrans", &twol, &blen, &negOnef, a->dG, &twol, a->theta, &one, &onef, a->sol, &one);
  RETURN (scs_int) info;
}

scs_int init_gels(Accel * a){
  DEBUG_FUNC
  blasint info;
  scs_float twork;
  scs_int l = a->l;
  blasint lwork = -1;
  blasint twol = 2 * l;
  blasint one = 1;
  blasint k = (blasint) a->k;
  blasint worksize;
  /* reuse variable names from posv */
  scs_free(a->X);
  scs_free(a->theta);
  a->theta = scs_malloc(sizeof(scs_float) * MAX(2 * a->l, a->k));
  a->X = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  BLAS(gels)("NoTrans", &twol, &k, &one, a->X, &twol, a->theta, &twol, &twork, &lwork, &info);
  worksize = (blasint) twork;
  a->worksize = worksize;
  a->wrk = scs_malloc(sizeof(blasint) * worksize);
  RETURN (scs_int) info;
}

scs_int solve_accel_linsys(Accel *a, scs_int len) {
  DEBUG_FUNC
  scs_int info;
  if (a->bad_numerics) {
    RETURN solve_with_gels(a, len);
  }
  info = solve_with_posv(a, len);
  if (info != 0) {
    scs_printf("\tposv ran into numerical issues, info %i, switching to gels\n", (int) info);
    a->bad_numerics = 1;
    init_gels(a);
    RETURN solve_with_gels(a, len);
  }
  RETURN info;
}

void update_accel_params(Work *w, scs_int idx) {
  DEBUG_FUNC
  scs_float *dF = w->accel->dF;
  scs_float *dG = w->accel->dG;
  scs_float *f = w->accel->f;
  scs_float *g = w->accel->g;
  scs_int l = w->m + w->n + 1;
  /* copy g_prev into idx col of dG */
  memcpy(&(dG[idx * 2 * l]), g, sizeof(scs_float) * 2 * l);
  /* copy f_prev into idx col of dF */
  memcpy(&(dF[idx * 2 * l]), f, sizeof(scs_float) * 2 * l);
  /* g = [u;v] */
  memcpy(g, w->u, sizeof(scs_float) * l);
  memcpy(&(g[l]), w->v, sizeof(scs_float) * l);
  /* calculcate f = g - [u_prev, v_prev] */
  memcpy(f, g, sizeof(scs_float) * 2 * l);
  addScaledArray(f, w->u_prev, l, -1.0);
  addScaledArray(&(f[l]), w->v_prev, l, -1.0);
  /* idx col of dG = g_prev - g */
  addScaledArray(&(dG[idx * 2 * l]), g, 2 * l, -1);
  /* idx col of dF = f_prev - f */
  addScaledArray(&(dF[idx * 2 * l]), f, 2 * l, -1);
  RETURN;
}

Accel *initAccel(Work *w) {
  DEBUG_FUNC
  Accel *a = scs_calloc(1, sizeof(Accel));
  if (!a) {
    RETURN SCS_NULL;
  }
  a->l = w->m + w->n + 1;
  a->k = MIN(2 * a->l, w->stgs->acceleration_lookback - 1);
  a->bad_numerics = 0;
  if (a->k <= 0) {
    RETURN a;
  }
  /* k = lookback - 1 since we use the difference form
     of anderson acceleration, and so there is one fewer var in lin sys. */
  /* Use MIN to prevent not full rank matrices */
  a->dF = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->dG = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->f = scs_calloc(2 * a->l, sizeof(scs_float));
  a->g = scs_calloc(2 * a->l, sizeof(scs_float));
  a->sol = scs_malloc(sizeof(scs_float) * 2 * a->l);
  /* posv swap to using these for gels if needed */
  a->theta = scs_malloc(sizeof(scs_float) * a->k);
  a->X = scs_malloc(sizeof(scs_float) * a->k * a->k);
  a->totalAccelTime = 0.0;
  if (!a->dF || !a->dG || !a->f || !a->g || !a->theta || !a->sol || !a->X) {
    freeAccel(a);
    a = SCS_NULL;
  }
  RETURN a;
}

scs_int accelerate(Work *w, scs_int iter) {
  DEBUG_FUNC
  scs_float *sol = w->accel->sol;
  scs_int l = w->accel->l;
  scs_int k = w->accel->k;
  scs_int info;
  timer accelTimer;
  tic(&accelTimer);
  if (k <= 0) {
    RETURN 0;
  }
  /* update dF, dG, f, g */
  update_accel_params(w, (k + iter - 1) % k);
  /* iter < k doesn't do any acceleration until iters hit k */
  /* if (iter < k) { */
  if (iter == 0) {
    RETURN 0;
  }
  /* solve linear system, new point stored in sol */
  info = solve_accel_linsys(w->accel, MIN(iter, k));
  /* set [u;v] = sol */
  if (info == 0) {
    memcpy(w->u, sol, sizeof(scs_float) * l);
    memcpy(w->v, &(sol[l]), sizeof(scs_float) * l);
  } else {
    /* Could fall back to gelsy or others in this case case */
    scs_printf("\taccelerate error, info %i, no acceleration applied\n", (int) info);
  }
  w->accel->totalAccelTime += tocq(&accelTimer);
  RETURN 0;
}

void freeAccel(Accel *a) {
  DEBUG_FUNC
  if (a) {
    if (a->dF) scs_free(a->dF);
    if (a->dG) scs_free(a->dG);
    if (a->f) scs_free(a->f);
    if (a->g) scs_free(a->g);
    if (a->sol) scs_free(a->sol);
    if (a->theta) scs_free(a->theta);
    if (a->X) scs_free(a->X);
    if (a->wrk) scs_free(a->wrk);
    scs_free(a);
  }
  RETURN;
}

#else

Accel *initAccel(Work *w) {
  Accel *a = scs_malloc(sizeof(Accel));
  a->totalAccelTime = 0.0;
  RETURN a;
}

void freeAccel(Accel *a) {
  if (a) {scs_free(a);}
}

scs_int accelerate(Work *w, scs_int iter) { RETURN 0; }
#endif

char *getAccelSummary(const Info *info, Accel *a) {
  DEBUG_FUNC
  char *str = scs_malloc(sizeof(char) * 64);
  sprintf(str, "\tAcceleration: avg solve time: %1.2es\n",
          a->totalAccelTime / (info->iter + 1) / 1e3);
  a->totalAccelTime = 0.0;
  RETURN str;
}
