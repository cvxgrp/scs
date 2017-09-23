#include "accel.h"
#include "scs.h"

void update_accel_params(Work * w, scs_int idx){
    DEBUG_FUNC
    scs_float* dF = w->accel->dF;
    scs_float* dG = w->accel->dG;
    scs_float* f = w->accel->f;
    scs_float* g = w->accel->g;
    scs_int l = w->m + w->n + 1;
    if (idx > 0) {
        // copy g_prev into idx col of dG
        memcpy(&(dG[(idx - 1) * 2 * l]), g, sizeof(scs_float) * 2 * l);
        // copy f_prev into idx col of dF
        memcpy(&(dF[(idx - 1) * 2 * l]), f, sizeof(scs_float) * 2 * l);
    }
    // g = [u;v]
    memcpy(g, w->u, sizeof(scs_float) * l);
    memcpy(&(g[l]), w->v, sizeof(scs_float) * l);
    // calulcate f = g - [u_prev, v_prev]
    memcpy(f, g, sizeof(scs_float) * 2 * l);
    addScaledArray(f, w->u_prev, l, -1.0);
    addScaledArray(&(f[l]), w->v_prev, l, -1.0);
    if (idx > 0) {
        // last col of dG = g_prev - g
        addScaledArray(&(dG[(idx - 1) * 2 * l]), g, 2 * l, -1);
        // last col of dF = f_prev - f
        addScaledArray(&(dF[(idx - 1) * 2 * l]), f, 2 * l, -1);
    }
    RETURN;
}

blasint initAccelWrk(Accel * a){
  blasint info;
  scs_float twork;
  scs_int l = a->l;
  blasint lwork = -1;
  blasint twol = 2 * l;
  blasint one = 1;
  scs_int k = a->k;
  blasint bk = (blasint) k;
  blasint worksize;
  BLAS(gels)("NoTrans", &twol, &bk, &one, a->dFQR, &twol, a->theta, &twol, &twork, &lwork, &info);
  worksize = (blasint) twork;
  a->wrk = scs_malloc(sizeof(blasint) * worksize);
  RETURN worksize;
}

Accel* initAccel(Work * w) {
  Accel * a = scs_malloc(sizeof(Accel));
  scs_int l = w->m + w->n + 1;
  a->l = l;
  /* k = lookback - 1 since we use the difference form
     of anderson acceleration, and so there is one fewer var in lin sys. */
  a->k = w->stgs->acceleration_lookback - 1;
  a->dF = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  a->dG = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  a->f = scs_malloc(sizeof(scs_float) * 2 * l);
  a->g = scs_malloc(sizeof(scs_float) * 2 * l);
  a->theta = scs_malloc(sizeof(scs_float) * 2 * l);
  a->tmp = scs_malloc(sizeof(scs_float) * 2 * l);
  a->dFQR = scs_malloc(sizeof(scs_float) * 2 * l * a->k);
  a->worksize = initAccelWrk(a);
  RETURN a;
}

void solve_accel_linsys(Accel * a) {
    blasint info;
    scs_int l = a->l;
    blasint twol = 2 * l;
    blasint one = 1;
    scs_int k = a->k;
    blasint bk = (blasint) k;
    scs_float negOnef = -1.0;
    scs_float onef = 1.0;
    memcpy(a->theta, a->f, sizeof(scs_float) * 2 * l);
    memcpy(a->dFQR, a->dF, sizeof(scs_float) * 2 * l * k);
    BLAS(gels)("NoTrans", &twol, &bk, &one, a->dFQR, &twol, a->theta, &twol, a->wrk, &(a->worksize), &info);
    memcpy(a->tmp, a->g, sizeof(scs_float) * 2 * l);
    // matmul g = g - dG * theta
    BLAS(gemv)("NoTrans", &twol, &bk, &negOnef, a->dG, &twol, a->theta, &one, &onef, a->tmp, &one);
}

scs_int accelerate(Work *w, scs_int iter) {
    // TODO: what if k > 2l ?
    DEBUG_FUNC
    scs_float* dF = w->accel->dF;
    scs_float* dG = w->accel->dG;
    scs_float * tmp = w->accel->tmp;
    scs_int l = w->accel->l;
    scs_int k = w->accel->k;
    if (k == 0) {
        RETURN 0;
    }
    if (iter < k) {
        update_accel_params(w, iter);
        RETURN 0;
    }
    if (iter > k) {
        // shift dF
        memcpy(dF, &(dF[2*l]), sizeof(scs_float) * 2 * l * (k-1));
        // shift dG
        memcpy(dG, &(dG[2*l]), sizeof(scs_float) * 2 * l * (k-1));
    }
    // update dF, dG, f, g
    update_accel_params(w, k);
    // solve linear system for theta
    solve_accel_linsys(w->accel);
    // set [u;v] = tmp
    memcpy(w->u, tmp, sizeof(scs_float) * l);
    memcpy(w->v, &(tmp[l]), sizeof(scs_float) * l);
    RETURN 0;
}

void freeAccel(Accel * a) {
    scs_free(a->dF);
    scs_free(a->dG);
    scs_free(a->f);
    scs_free(a->g);
    scs_free(a->dFQR);
    scs_free(a->wrk);
    scs_free(a->theta);
    scs_free(a->tmp);
    scs_free(a);
}

