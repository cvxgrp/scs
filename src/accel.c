#include "accel.h"
#include "scs.h"
#include "scs_blas.h"

struct SCS_ACCEL {
#ifdef LAPACK_LIB_FOUND
    scs_float * dF;
    scs_float * dG;
    scs_float * f;
    scs_float * g;
    scs_float * theta;
    scs_float * tmp;
    scs_float * dFQR;
    scs_float * wrk;
    scs_int k, l;
    blasint worksize;
#endif
    scs_float totalAccelTime;
};

#ifdef LAPACK_LIB_FOUND
void BLAS(gels)(const char *trans, const blasint *m, const blasint *n, const blasint *nrhs, scs_float *a, const blasint *lda, scs_float *b, const blasint *ldb, scs_float *work, const blasint *lwork, blasint *info);
void BLAS(gemv)(const char *trans, const blasint *m, const blasint *n, const scs_float *alpha, const scs_float *a, const blasint *lda, const scs_float *x, const blasint *incx, const scs_float *beta, scs_float *y, const blasint *incy);

blasint initAccelWrk(Accel * a);
scs_int solve_accel_linsys(Accel * a);
void update_accel_params(Work * w, scs_int idx);

void update_accel_params(Work * w, scs_int idx){
    DEBUG_FUNC
    scs_float* dF = w->accel->dF;
    scs_float* dG = w->accel->dG;
    scs_float* f = w->accel->f;
    scs_float* g = w->accel->g;
    scs_int l = w->m + w->n + 1;
    if (idx > 0) {
        /* copy g_prev into idx col of dG */
        memcpy(&(dG[(idx - 1) * 2 * l]), g, sizeof(scs_float) * 2 * l);
        /* copy f_prev into idx col of dF */
        memcpy(&(dF[(idx - 1) * 2 * l]), f, sizeof(scs_float) * 2 * l);
    }
    /* g = [u;v] */
    memcpy(g, w->u, sizeof(scs_float) * l);
    memcpy(&(g[l]), w->v, sizeof(scs_float) * l);
    /* calulcate f = g - [u_prev, v_prev] */
    memcpy(f, g, sizeof(scs_float) * 2 * l);
    addScaledArray(f, w->u_prev, l, -1.0);
    addScaledArray(&(f[l]), w->v_prev, l, -1.0);
    if (idx > 0) {
        /* last col of dG = g_prev - g */
        addScaledArray(&(dG[(idx - 1) * 2 * l]), g, 2 * l, -1);
        /* last col of dF = f_prev - f */
        addScaledArray(&(dF[(idx - 1) * 2 * l]), f, 2 * l, -1);
    }
    RETURN;
}

blasint initAccelWrk(Accel * a){
  DEBUG_FUNC
  blasint info;
  scs_float twork;
  scs_int l = a->l;
  blasint lwork = -1;
  blasint twol = 2 * l;
  blasint one = 1;
  blasint k = (blasint) a->k;
  blasint worksize;
  BLAS(gels)("NoTrans", &twol, &k, &one, a->dFQR, &twol, a->theta, &twol, &twork, &lwork, &info);
  worksize = (blasint) twork;
  a->worksize = worksize;
  a->wrk = scs_malloc(sizeof(blasint) * worksize);
  RETURN info;
}

Accel* initAccel(Work * w) {
  DEBUG_FUNC
  Accel * a = scs_malloc(sizeof(Accel));
  scs_int l = w->m + w->n + 1;
  scs_int info;
  if (!a) {
    RETURN SCS_NULL;
  }
  a->l = l;
  /* k = lookback - 1 since we use the difference form
     of anderson acceleration, and so there is one fewer var in lin sys. */
  /* Use MIN to prevent not full rank matrices */
  a->k = MIN(2 * l, w->stgs->acceleration_lookback - 1);
  if (a->k <= 0) {
    a->dF = a->dG = a->f = a->g = a->theta = a->tmp = a->dFQR = a->wrk = SCS_NULL;
    RETURN a;
  }
  a->dF = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  a->dG = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  a->f = scs_malloc(sizeof(scs_float) * 2 * l);
  a->g = scs_malloc(sizeof(scs_float) * 2 * l);
  a->theta = scs_malloc(sizeof(scs_float) * MAX(2 * l, a->k));
  a->tmp = scs_malloc(sizeof(scs_float) * 2 * l);
  a->dFQR = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  a->totalAccelTime = 0.0;
  info = initAccelWrk(a);
  if (!a->dF || ! a->dG || !a->f || !a->g || !a->theta || !a->tmp || ! a->dFQR || info != 0) {
    freeAccel(a);
    a = SCS_NULL;
  }
  RETURN a;
}

scs_int solve_accel_linsys(Accel * a) {
  DEBUG_FUNC
  blasint info;
  scs_int l = a->l;
  blasint twol = 2 * l;
  blasint one = 1;
  blasint k = (blasint) a->k;
  scs_float negOnef = -1.0;
  scs_float onef = 1.0;
  memcpy(a->theta, a->f, sizeof(scs_float) * 2 * l);
  memcpy(a->dFQR, a->dF, sizeof(scs_float) * 2 * l * k);
  BLAS(gels)("NoTrans", &twol, &k, &one, a->dFQR, &twol, a->theta, &twol, a->wrk, &(a->worksize), &info);
  memcpy(a->tmp, a->g, sizeof(scs_float) * 2 * l);
  /* matmul g = g - dG * theta */
  BLAS(gemv)("NoTrans", &twol, &k, &negOnef, a->dG, &twol, a->theta, &one, &onef, a->tmp, &one);
  RETURN info;
}

scs_int accelerate(Work *w, scs_int iter) {
  DEBUG_FUNC
  scs_float* dF = w->accel->dF;
  scs_float* dG = w->accel->dG;
  scs_float * tmp = w->accel->tmp;
  scs_int l = w->accel->l;
  scs_int k = w->accel->k;
  scs_int info;
  timer accelTimer;
  tic(&accelTimer);
  if (k <= 0) {
    RETURN 0;
  }
  if (iter < k) {
    update_accel_params(w, iter);
    RETURN 0;
  }
  if (iter > k) {
    /* shift dF */
    memmove(dF, &(dF[2*l]), sizeof(scs_float) * 2 * l * (k-1));
    /* shift dG */
    memmove(dG, &(dG[2*l]), sizeof(scs_float) * 2 * l * (k-1));
  }
  /* update dF, dG, f, g */
  update_accel_params(w, k);
  /* solve linear system for theta */
  info = solve_accel_linsys(w->accel);
  /* set [u;v] = tmp */
  memcpy(w->u, tmp, sizeof(scs_float) * l);
  memcpy(w->v, &(tmp[l]), sizeof(scs_float) * l);
  w->accel->totalAccelTime += tocq(&accelTimer);
  RETURN info;
}

void freeAccel(Accel * a) {
  DEBUG_FUNC
    if(a) {
      if (a->dF)
        scs_free(a->dF);
      if (a->dG)
        scs_free(a->dG);
      if (a->f)
        scs_free(a->f);
      if (a->g)
        scs_free(a->g);
      if (a->dFQR)
        scs_free(a->dFQR);
      if (a->theta)
        scs_free(a->theta);
      if (a->tmp)
        scs_free(a->tmp);
      if (a->wrk)
        scs_free(a->wrk);
      scs_free(a);
    }
  RETURN;
}

#else

Accel* initAccel(Work * w) {
  Accel * a = scs_malloc(sizeof(Accel));
  a->totalAccelTime = 0.0;
  RETURN a;
}

void freeAccel(Accel * a) {
  if (a)
    scs_free(a);
}

scs_int accelerate(Work *w, scs_int iter) {
 RETURN 0;
}
#endif

char *getAccelSummary(const Info *info, Accel *a) {
  DEBUG_FUNC
    char *str = scs_malloc(sizeof(char) * 64);
  sprintf(str, "\tAcceleration: avg solve time: %1.2es\n",
      a->totalAccelTime / (info->iter + 1) / 1e3);
  a->totalAccelTime = 0.0;
  RETURN str;
}
