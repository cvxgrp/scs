#ifndef ACCEL_H_GUARD
#define ACCEL_H_GUARD

#include "glbopts.h"
#include "scs_blas.h"

void BLAS(gels)(const char *trans, const blasint *m, const blasint *n, const blasint *nrhs, scs_float *a, const blasint *lda, scs_float *b, const blasint *ldb, scs_float *work, const blasint *lwork, blasint *info);
void BLAS(gemv)(const char *trans, const blasint *m, const blasint *n, const scs_float *alpha, const scs_float *a, const blasint *lda, const scs_float *x, const blasint *incx, const scs_float *beta, scs_float *y, const blasint *incy);

typedef struct {
    scs_float * dF;
    scs_float * dG;
    scs_float * f;
    scs_float * g;
    scs_float * theta;
    scs_float * tmp;
    scs_int k, l;
    scs_float * dFQR;
    scs_float * wrk;
    blasint worksize;
} Accel;

void update_accel_params(Work * w, scs_int idx);
Accel* initAccel(Work * w);
blasint initAccelWrk(Accel * a);
void freeAccel(Accel * a);
scs_int accelerate(Work *w, scs_int iter);
void solve_accel_linsys(Accel * a);

#endif
