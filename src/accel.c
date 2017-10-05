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
  scs_float *scratch;
  //scs_float *Q;
  scs_float *dFold;
  scs_float *R;
  scs_float *delta;
  scs_int k, l;
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
void BLAS(geqrf)(blasint *m, blasint *n, scs_float *a, blasint * lda, scs_float *tau, scs_float *work, blasint *lwork, blasint *info);
void BLAS(orgqr)(blasint *m, blasint *n, blasint *k, scs_float * a, blasint *lda, scs_float *tau, scs_float *work, blasint *lwork, blasint *info);
void BLAS(trsv)(const char *uplo, const char *trans, const char *diag,
                blasint *n, scs_float *a, blasint *lda, scs_float *x,
                blasint *incx);
scs_float BLAS(nrm2)(const blasint *n, scs_float *x, const blasint *incx);
void BLAS(rotg) (scs_float *a, scs_float *b, scs_float *c, scs_float *s);
void BLAS(rot) (const blasint *n, scs_float *x, const blasint *incx, scs_float *y, const blasint *incy, const scs_float *c, const scs_float *s);

void BLAS(gels)(const char *trans, const blasint *m, const blasint *n,
                const blasint *nrhs, scs_float *a, const blasint *lda,
                scs_float *b, const blasint *ldb, scs_float *work,
                const blasint *lwork, blasint *info);
void BLAS(gemm)(const char *transa, const char *transb, blasint *m, blasint *
    n, blasint *k, scs_float *alpha, scs_float *a, blasint *lda,
    scs_float *b, blasint *ldb, scs_float *beta, scs_float *c, blasint
    *ldc);

/* dF * sol = QR * sol = f */
scs_float * solve_accel_linsys(Accel *a, scs_int len) {
  DEBUG_FUNC
  blasint twol = 2 * a->l;
  blasint one = 1;
  blasint bk = (blasint) a->k;
  scs_float onef = 1.0;
  scs_float zerof = 0.0;
  scs_float negOnef = -1.0;
  /* sol = f */
  memcpy(a->sol, a->f, sizeof(scs_float) * 2 * a->l);
  /* scratch = Q' * sol = R^-T * dF' * sol*/
  BLAS(gemv)("Trans", &twol, &bk, &onef, a->dF, &twol, a->sol, &one, &zerof, a->scratch, &one);
  BLAS(trsv)("Upper", "Trans", "NotUnitDiag", &bk, a->R, &bk, a->scratch, &one);
  /* scratch = R^-1 * scratch */
  BLAS(trsv)("Upper", "NoTrans", "NotUnitDiag", &bk, a->R, &bk, a->scratch, &one);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * a->l);
  /* sol = sol - dG * scratch */
  BLAS(gemv)("NoTrans", &twol, &bk, &negOnef, a->dG, &twol, a->scratch, &one, &onef, a->sol, &one);
  RETURN a->sol;
}

scs_int solve_with_gels(Accel * a, scs_int len) {
  DEBUG_FUNC
  blasint info;
  scs_int l = a->l;
  blasint twol = 2 * l;
  blasint one = 1;
  blasint blen = (blasint) len;
  scs_float negOnef = -1.0;
  scs_float onef = 1.0;
  scs_float * X = scs_malloc(sizeof(scs_float) * 2 * l * len);
  memcpy(a->scratch, a->f, sizeof(scs_float) * 2 * l);
  memcpy(X, a->dF, sizeof(scs_float) * 2 * l * len);
  scs_float * wrk = scs_malloc(sizeof(scs_float) * 2 * l * len);
  blasint worksize = 2 * l * len;
  /* solve dF scratch = f */
  BLAS(gels)("NoTrans", &twol, &blen, &one, X, &twol, a->scratch, &twol, wrk, &worksize, &info);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * l);
  /* sol = sol - dG * scratch */
  BLAS(gemv)("NoTrans", &twol, &blen, &negOnef, a->dG, &twol, a->scratch, &one, &onef, a->sol, &one);
  scs_free(wrk);
  scs_free(X);
  RETURN info;
}

void update_accel_params(Work *w, scs_int idx) {
  DEBUG_FUNC
  scs_float *dF = w->accel->dF;
  scs_float *dG = w->accel->dG;
  scs_float *f = w->accel->f;
  scs_float *g = w->accel->g;
  scs_float *delta = w->accel->delta;
  Accel *a = w->accel;
  scs_int l = w->m + w->n + 1;
  memcpy(a->dFold, dF, sizeof(scs_float) * 2 * l * a->k);
  /* copy old col into delta */
  memcpy(delta, &(dF[idx * 2 * l]), sizeof(scs_float) * 2 * l);
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
  /* delta = delta - new col */
  addScaledArray(delta, &(dF[idx * 2 * l]), 2 * l, -1.0);
  /* delta = new - old */
  scaleArray(delta, -1.0, 2 * l);
  //scs_printf("norm f %e\n", calcNorm(f, 2 * l));
  //scs_printf("norm delta %e\n", calcNorm(delta, 2 * l));
  RETURN;
}

Accel *initAccel(Work *w) {
  DEBUG_FUNC
  Accel *a = scs_calloc(1, sizeof(Accel));
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
  a->dF = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->dFold = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->dG = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->f = scs_calloc(2 * a->l, sizeof(scs_float));
  a->g = scs_calloc(2 * a->l, sizeof(scs_float));
  //a->Q = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->R = scs_calloc(a->k * a->k, sizeof(scs_float));
  a->sol = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->scratch = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->delta = scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->totalAccelTime = 0.0;
  if (!a->dF || !a->dG || !a->f || !a->g || !a->scratch || !a->sol ||
      !a->delta || !a->R) {
    freeAccel(a);
    a = SCS_NULL;
  }
  RETURN a;
}

void qrfactorize(Accel * a) {
  DEBUG_FUNC
  scs_int l = a->l;
  scs_int i;
  blasint twol = 2 * l;
  blasint bk = (blasint) a->k;
  blasint negOne = -1;
  blasint info;
  blasint lwork;
  scs_float worksize;
  scs_float * work;
  scs_float * tau = scs_malloc(a->k * sizeof(scs_float));
  scs_float * Q = scs_malloc(2 * a->l * a->k * sizeof(scs_float));
  memcpy(Q, a->dF, sizeof(scs_float) * a->k * 2 * a->l);
  BLAS(geqrf)(&twol, &bk, Q, &twol, tau, &worksize, &negOne, &info);
  lwork = (blasint) worksize;
  work = scs_malloc(lwork * sizeof(scs_float));
  BLAS(geqrf)(&twol, &bk, Q, &twol, tau, work, &lwork, &info);
  //scs_printf("info %i\n", info);
  scs_free(work);
  //scs_printf("norm dF %e\n", calcNorm(a->dF, 2 * a->l * a->k));
  for (i = 0; i < a->k; ++i) {
    memcpy(&(a->R[i * a->k]), &(Q[i * a->l * 2]), sizeof(scs_float) * (i + 1));
  }
  //BLAS(orgqr)(&twol, &bk, &bk, a->Q, &twol, tau, &worksize, &negOne, &info);
  //lwork = (blasint) worksize;
  //work = scs_malloc(lwork * sizeof(scs_float));
  //BLAS(orgqr)(&twol, &bk, &bk, a->Q, &twol, tau, work, &lwork, &info);
  //scs_free(work);
  //for (i = 0; i < a->k; ++i){
  //  scs_printf("R[%i, %i] = %e, ", i, i, a->R[i * a->k + i]);
  //}
  scs_free(Q);
  scs_free(tau);
  RETURN;
}

void update_factorization(Accel * a, scs_int idx) {
  DEBUG_FUNC
  //scs_float * Q = a->Q;
  scs_float * R = a->R;
  scs_float * u = a->delta;
  scs_float * w = a->scratch;
  blasint one = 1;
  blasint bk = (blasint) a->k;
  blasint twol = (blasint) 2 * a->l;
  scs_float zerof = 0.0;
  scs_float onef = 1.0;
  scs_float negOnef = -1.0;
  scs_float nrm_u;
  scs_float c;
  scs_float s;

  scs_int k = a->k;
  scs_int l = a->l;
  scs_int i;

  /* w = Q' * delta = R^-T * dF' * delta, size k: col of R */
  BLAS(gemv)("Trans", &twol, &bk, &onef, a->dFold, &twol, u, &one, &zerof, w, &one);
  BLAS(trsv)("Upper", "Trans", "NotUnitDiag", &bk, a->R, &bk, w, &one);
  //BLAS(gemv)("Trans", &twol, &bk, &onef, a->Q, &twol, u, &one, &zerof, w, &one);
  /* u = delta - Q * w = dF * R^-1 w, size m: col of Q */
  scs_float * tmp = scs_malloc(sizeof(scs_float) * k);
  memcpy(tmp, w, sizeof(scs_float) * k);
  BLAS(trsv)("Upper", "NoTrans", "NotUnitDiag", &bk, a->R, &bk, tmp, &one);
  BLAS(gemv)("NoTrans", &twol, &bk, &negOnef, a->dFold, &twol, tmp, &one, &onef, u, &one);
  scs_free(tmp);
  /* nrm_u = ||u|| */
  nrm_u = BLAS(nrm2)(&twol, u, &one);
  /* u = u / ||u|| */
  scaleArray(u, 1.0 / nrm_u, 2 * l);
  /* R col += w */
  addScaledArray(&(R[idx * k]), w, k, 1.0);

  scs_float * bot_row = scs_calloc(k, sizeof(scs_float));
  /* Givens rotations, start with fake bottom row of R, extra col of Q */
  //scs_printf("idx %i\n", idx);
  //scs_printf("k %i\n", k);
  scs_int ridx = k * idx + k - 1;
  scs_float r1 = R[ridx];
  scs_float r2 = nrm_u;
  bot_row[idx] = nrm_u;
  BLAS(rotg)(&r1, &r2, &c, &s);
  BLAS(rot)(&bk, &(R[k - 1]), &bk, bot_row, &one, &c, &s);
  //BLAS(rot)(&twol, &(Q[2 * l * (k - 1)]), &one, u, &one, &c, &s);

  /* Walk up the spike, R finishes upper Hessenberg */
  for (i = k; i > idx + 1; --i) {
    scs_int ridx = k * idx + i - 1;
    scs_float r1 = R[ridx - 1];
    scs_float r2 = R[ridx];
    BLAS(rotg)(&(r1), &(r2), &c, &s);
    BLAS(rot)(&bk, &(R[i - 2]), &bk, &(R[i - 1]), &bk, &c, &s);
    //BLAS(rot)(&twol, &(Q[2 * l * (i-2)]), &one, &(Q[2 * l * (i-1)]), &one, &c, &s);
  }

  /* Walk down the sub-diagonal, R finishes upper triangular */
  for (i = idx + 1; i < k - 1; ++i) {
    scs_int ridx = k * i + i;
    scs_float r1 = R[ridx];
    scs_float r2 = R[ridx + 1];
    BLAS(rotg)(&r1, &r2, &c, &s);
    BLAS(rot)(&bk, &(R[i]), &bk, &(R[i + 1]), &bk, &c, &s);
    //BLAS(rot)(&twol, &(Q[2 * l * i]), &one, &(Q[2 * l * (i+1)]), &one, &c, &s);
  }

  scs_float min_r = 9999999.9;
  scs_float prod_r = 1.0;
  for (i = 0; i < k; ++i){
    //scs_printf("R[%i, %i] = %e, ", i, i, R[i * k + i]);
    /*
    if (ABS(R[i * k + i]) < min_r){
        min_r = ABS(R[i * k + i]);
        prod_r *= ABS(R[i * k + i]);
    }

    if (R[i * k + i] < 0) {
      R[i * k + i] -= 1e-1;
    } else {
      R[i * k + i] += 1e-1;
    }
    */
  }
  //scs_printf("\n");
  //scs_printf("min diag R %e\n", min_r);
  //scs_printf("prod diag R %e\n", prod_r);

  /* Finish fake bottom row of R, extra col of Q */
  BLAS(rotg)(&(R[k * k - 1]), &(bot_row[k - 1]), &c, &s);
  //BLAS(rot)(&twol, &(Q[2 * l * (k - 1)]), &one, u, &one, &c, &s);
  scs_free(bot_row);
  RETURN;
}

scs_int accelerate(Work *w, scs_int iter) {
  DEBUG_FUNC
  //scs_printf("iter %i\n", iter);
  scs_int l = w->accel->l;
  scs_int k = w->accel->k;
  scs_float *sol;
  Accel * a = w->accel;
  timer accelTimer;
  tic(&accelTimer);
  if (k <= 0) {
    RETURN 0;
  }
  /* update dF, dG, f, g, delta */
  update_accel_params(w, (k + iter - 1) % k);
  /* iter < k doesn't do any acceleration until iters hit k */
  if (iter < k) {
    RETURN 0;
  }
  if (iter == k) {
    qrfactorize(w->accel);
    //printArray(w->accel->Q, k * 2 * l, "Q_true");
    //printArray(w->accel->R, k * k, "R_true");
  } else {
    /* update Q, R factors */

    update_factorization(w->accel, (k + iter - 1) % k);
    //qrfactorize(w->accel);
    //scs_float * dF0 = scs_calloc(2 * a->l * a->k, sizeof(scs_float));
    //blasint twol = 2 * a->l;
    //blasint bk = (blasint) a->k;
    //scs_float onef = 1.0;
    //scs_float zerof = 0.0;


    //BLAS(gemm)("NoTrans", "NoTrans", &twol, &bk, &bk, &onef, a->Q, &twol,
    //    a->R, &bk, &zerof, dF0, &twol);
    //printArray(w->accel->Q, k * 2 * l, "Q");
    //printArray(w->accel->R, k * k, "R");

    //scs_printf("||DdF|| = %e\n", calcNormDiff(a->dF, dF0, 2 * a->l * a->k));
    //scs_printf("||dF|| = %e\n", calcNorm(a->dF, 2 * a->l * a->k));
    //scs_printf("||DdF||/||dF|| = %e\n", calcNormDiff(a->dF, dF0, 2 * a->l *
    //      a->k) / calcNorm(a->dF, 2 * a->l * a->k));
    //scs_free(dF0);
    //printArray(w->accel->Q, k * 2 * l, "Q_true");
    //printArray(w->accel->R, k * k, "R_true");

  }
  /* solve linear system, new point stored in sol */
  sol = solve_accel_linsys(w->accel, k);
  //printArray(sol, 2 * l, "sol1");
  //solve_with_gels(w->accel, k);
  //printArray(sol, 2 * l, "sol2");
  /* set [u;v] = sol */
  //scs_printf("norm of sol %e\n", calcNorm(sol, 2 *l));
  memcpy(w->u, sol, sizeof(scs_float) * l);
  memcpy(w->v, &(sol[l]), sizeof(scs_float) * l);
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
    if (a->scratch) scs_free(a->scratch);
    if (a->dFold) scs_free(a->dFold);
    //if (a->Q) scs_free(a->Q);
    if (a->R) scs_free(a->R);
    if (a->delta) scs_free(a->delta);
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
