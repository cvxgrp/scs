/*
 * Anderson acceleration.
 *
 * x: input iterate
 * x_prev: previous input iterate
 * f: f(x) output of map f applied to x
 * g: x - f error
 * g_prev: previous error
 * s: x - x_prev
 * y: g - g_prev
 * d: s - y
 *
 * capital letters are the variables stacked columnwise
 * idx tracks current index where latest quantities written
 * idx cycles from left to right columns in matrix
 *
 * Type-I:
 * return x - (S - Y) * ( S'Y + r I)^{-1} ( S'g )
 *
 * Type-II:
 * return x - (S - Y) * ( Y'Y + r I)^{-1} ( Y'g )
 *
 * TODO: need to fully implement safeguards
 *
 */

#include "aa.h"
#include "scs_blas.h"

#if PROFILING > 0

#define TIME_TIC \
  timer __t;     \
  tic(&__t);
#define TIME_TOC toc(__func__, &__t);

#include <time.h>
typedef struct timer {
  struct timespec tic;
  struct timespec toc;
} timer;

void tic(timer *t) { clock_gettime(CLOCK_MONOTONIC, &t->tic); }

aa_float tocq(timer *t) {
  struct timespec temp;

  clock_gettime(CLOCK_MONOTONIC, &t->toc);

  if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
    temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
  } else {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
    temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
  }
  return (aa_float)temp.tv_sec * 1e3 + (aa_float)temp.tv_nsec / 1e6;
}

aa_float toc(const char *str, timer *t) {
  aa_float time = tocq(t);
  scs_printf("%s - time: %8.4f milli-seconds.\n", str, time);
  return time;
}

#else

#define TIME_TIC
#define TIME_TOC

#endif


/* BLAS functions used */
aa_float BLAS(nrm2)(blas_int *n, aa_float *x, blas_int *incx);
void BLAS(axpy)(blas_int *n, aa_float *a, const aa_float *x, blas_int *incx,
                aa_float *y, blas_int *incy);
void BLAS(gemv)(const char *trans, const blas_int *m, const blas_int *n,
                const aa_float *alpha, const aa_float *a, const blas_int *lda,
                const aa_float *x, const blas_int *incx, const aa_float *beta,
                aa_float *y, const blas_int *incy);
void BLAS(gesv)(blas_int *n, blas_int *nrhs, aa_float *a, blas_int *lda,
                blas_int *ipiv, aa_float *b, blas_int *ldb, blas_int *info);
void BLAS(gemm)(const char *transa, const char *transb, blas_int *m,
                blas_int *n, blas_int *k, aa_float *alpha, aa_float *a,
                blas_int *lda, aa_float *b, blas_int *ldb, aa_float *beta,
                aa_float *c, blas_int *ldc);


/* This file uses Anderson acceleration to improve the convergence of
 * a fixed point mapping.
 * At each iteration we need to solve a (small) linear system, we
 * do this using LAPACK ?gesv.
 */

/* contains the necessary parameters to perform aa at each step */
struct ACCEL_WORK {
  aa_int type1; /* bool, if true type 1 aa otherwise type 2 */
  aa_int mem;   /* aa memory */
  aa_int dim;   /* variable dimension */
  aa_int iter;  /* current iteration */

  aa_float eta; /* regularization */

  aa_float *x; /* x input to map*/
  aa_float *f; /* f(x) output of map */
  aa_float *g; /* x - f(x) */

  /* from previous iteration */
  aa_float *g_prev; /* x - f(x) */

  aa_float *y; /* g - g_prev */
  aa_float *s; /* x - x_prev */
  aa_float *d; /* f - f_prev */

  aa_float *Y; /* matrix of stacked y values */
  aa_float *S; /* matrix of stacked s values */
  aa_float *D; /* matrix of stacked d values = (S-Y) */
  aa_float *M; /* S'Y or Y'Y depending on type of aa */

  /* workspace variables */
  aa_float *work; /* scratch space */
  blas_int *ipiv; /* permutation variable, not used after solve */
};

/* sets a->M to S'Y or Y'Y depending on type of aa used */
/* M is len x len after this */
static void set_m(AaWork *a, aa_int len) {
  TIME_TIC
  aa_float r, nrm_y, nrm_s; /* add r to diags for regularization */
  aa_int i;
  blas_int bdim = (blas_int)(a->dim), one = 1;
  blas_int blen = (blas_int)len, btotal = (blas_int)(a->dim * len);
  aa_float onef = 1.0, zerof = 0.0;
  /* if len < mem this only uses len cols */
  BLAS(gemm)("Trans", "No", &blen, &blen, &bdim, &onef, a->type1 ? a->S : a->Y,
              &bdim, a->Y, &bdim, &zerof, a->M, &blen);
  if (a->eta > 0) {
    /* TODO: this regularization doesn't make much sense for type-I */
    /* but we do it anyway since it seems to help */
    /* typically type-I does better with higher eta (more reg) than type-II */
    nrm_y = BLAS(nrm2)(&btotal, a->Y, &one);
    nrm_s = BLAS(nrm2)(&btotal, a->S, &one);
    r = a->eta * (nrm_y * nrm_y + nrm_s * nrm_s);
    #if EXTRA_VERBOSE > 5
    scs_printf("iter: %i, len: %i, norm: Y %.2e, norm: S %.2e, r: %.2e\n", a->iter,
            len, nrm_y, nrm_s, r);
    #endif
    for (i = 0; i < len; ++i){
      a->M[i + len * i] += r;
    }
  }
  TIME_TOC
  return;
}

/* initialize accel params, in particular x_prev, f_prev, g_prev */
static void init_accel_params(const aa_float *x, const aa_float *f,
                              AaWork *a) {
  blas_int bdim = (blas_int)a->dim;
  aa_float neg_onef = -1.0;
  blas_int one = 1;
  /* x_prev = x */
  memcpy(a->x, x, sizeof(aa_float) * a->dim);
  /* f_prev = f */
  memcpy(a->f, f, sizeof(aa_float) * a->dim);
  /* g_prev = x */
  memcpy(a->g_prev, x, sizeof(aa_float) * a->dim);
  /* g_prev = x_prev - f_prev */
  BLAS(axpy)(&bdim, &neg_onef, f, &one, a->g_prev, &one);
}

/* updates the workspace parameters for aa for this iteration */
static void update_accel_params(const aa_float *x, const aa_float *f,
                                AaWork *a, aa_int len) {
  /* at the start a->x = x_prev and a->f = f_prev */
  TIME_TIC
  aa_int idx = (a->iter - 1) % a->mem;
  blas_int one = 1;
  blas_int bdim = (blas_int)a->dim;
  aa_float neg_onef = -1.0;

  /* g = x */
  memcpy(a->g, x, sizeof(aa_float) * a->dim);
  /* s = x */
  memcpy(a->s, x, sizeof(aa_float) * a->dim);
  /* d = f */
  memcpy(a->d, f, sizeof(aa_float) * a->dim);
  /* g =  x - f */
  BLAS(axpy)(&bdim, &neg_onef, f, &one, a->g, &one);
  /* s = x - x_prev */
  BLAS(axpy)(&bdim, &neg_onef, a->x, &one, a->s, &one);
  /* d = f - f_prev */
  BLAS(axpy)(&bdim, &neg_onef, a->f, &one, a->d, &one);

  /* g, s, d correct here */

  /* y = g */
  memcpy(a->y, a->g, sizeof(aa_float) * a->dim);
  /* y = g - g_prev */
  BLAS(axpy)(&bdim, &neg_onef, a->g_prev, &one, a->y, &one);

  /* y correct here */

  /* copy y into idx col of Y */
  memcpy(&(a->Y[idx * a->dim]), a->y, sizeof(aa_float) * a->dim);
  /* copy s into idx col of S */
  memcpy(&(a->S[idx * a->dim]), a->s, sizeof(aa_float) * a->dim);
  /* copy d into idx col of D */
  memcpy(&(a->D[idx * a->dim]), a->d, sizeof(aa_float) * a->dim);

  /* Y, S, D correct here */

  memcpy(a->f, f, sizeof(aa_float) * a->dim);
  memcpy(a->x, x, sizeof(aa_float) * a->dim);

  /* x, f correct here */

  /* set M = S'Y or Y'Y depending on type of aa used */
  set_m(a, len);

  /* M correct here */

  memcpy(a->g_prev, a->g, sizeof(aa_float) * a->dim);

  /* g_prev set for next iter here */

  TIME_TOC
  return;
}

/* solves the system of equations to perform the aa update
 * at the end f contains the next iterate to be returned
 */
static aa_int solve(aa_float *f, AaWork *a, aa_int len) {
  TIME_TIC
  blas_int info = -1, bdim = (blas_int)(a->dim), one = 1, blen = (blas_int)len;
  blas_int bmem = (blas_int)a->mem;
  aa_float neg_onef = -1.0, onef = 1.0, zerof = 0.0, nrm;
  /* work = S'g or Y'g */
  BLAS(gemv)("Trans", &bdim, &blen, &onef, a->type1 ? a->S : a->Y, &bdim, a->g,
              &one, &zerof, a->work, &one);
  /* work = M \ work, where update_accel_params has set M = S'Y or M = Y'Y */
  BLAS(gesv)(&blen, &one, a->M, &blen, a->ipiv, a->work, &blen, &info);
  nrm = BLAS(nrm2)(&bmem, a->work, &one);
  #if EXTRA_VERBOSE > 10
  scs_printf("AA type %i, iter: %i, len %i, info: %i, norm %.2e\n",
          a->type1 ? 1 : 2, (int)a->iter, (int) len, (int)info, nrm);
  #endif
  if (info < 0 || nrm >= MAX_AA_NRM) {
    #if EXTRA_VERBOSE > 0
    scs_printf("Error in AA type %i, iter: %i, len %i, info: %i, norm %.2e\n",
            a->type1 ? 1 : 2, (int)a->iter, (int) len, (int)info, nrm);
    #endif
    /* reset aa for stability */
    aa_reset(a);
    TIME_TOC
    return -1;
  }
  /* if solve was successful then set f -= D * work */
  BLAS(gemv)("NoTrans", &bdim, &blen, &neg_onef, a->D, &bdim, a->work, &one,
             &onef, f, &one);
  TIME_TOC
  return (aa_int)info;
}

/*
 * API functions below this line, see aa.h for descriptions.
 */
AaWork *aa_init(aa_int dim, aa_int mem, aa_int type1, aa_float eta) {
  AaWork *a = (AaWork *)calloc(1, sizeof(AaWork));
  if (!a) {
    scs_printf("Failed to allocate memory for AA.\n");
    return (void *)0;
  }
  a->type1 = type1;
  a->iter = 0;
  a->dim = dim;
  a->mem = mem;
  a->eta = eta;
  if (a->mem <= 0) {
    return a;
  }

  a->x = (aa_float *)calloc(a->dim, sizeof(aa_float));
  a->f = (aa_float *)calloc(a->dim, sizeof(aa_float));
  a->g = (aa_float *)calloc(a->dim, sizeof(aa_float));

  a->g_prev = (aa_float *)calloc(a->dim, sizeof(aa_float));

  a->y = (aa_float *)calloc(a->dim, sizeof(aa_float));
  a->s = (aa_float *)calloc(a->dim, sizeof(aa_float));
  a->d = (aa_float *)calloc(a->dim, sizeof(aa_float));

  a->Y = (aa_float *)calloc(a->dim * a->mem, sizeof(aa_float));
  a->S = (aa_float *)calloc(a->dim * a->mem, sizeof(aa_float));
  a->D = (aa_float *)calloc(a->dim * a->mem, sizeof(aa_float));

  a->M = (aa_float *)calloc(a->mem * a->mem, sizeof(aa_float));
  a->work = (aa_float *)calloc(a->mem, sizeof(aa_float));
  a->ipiv = (blas_int *)calloc(a->mem, sizeof(blas_int));
  return a;
}

aa_int aa_apply(aa_float *f, const aa_float *x, AaWork *a) {
  TIME_TIC
  aa_int status;
  aa_int len = MIN(a->iter, a->mem);
  if (a->mem <= 0) {
    return 0;
  }
  if (a->iter == 0) {
    /* if first iteration then seed params for next iter */
    init_accel_params(x, f, a);
    a->iter++;
    TIME_TOC
    return 0;
  }
  /* set various accel quantities */
  update_accel_params(x, f, a, len);
  /* solve linear system, new point overwrites f if successful */
  status = solve(f, a, len);
  a->iter++;
  TIME_TOC
  return status;
}

void aa_finish(AaWork *a) {
  if (a) {
    free(a->x);
    free(a->f);
    free(a->g);
    free(a->g_prev);
    free(a->y);
    free(a->s);
    free(a->d);
    free(a->Y);
    free(a->S);
    free(a->D);
    free(a->M);
    free(a->work);
    free(a->ipiv);
    free(a);
  }
  return;
}

void aa_reset(AaWork *a) {
  /* to reset we simply set a->iter = 0 */
  a->iter = 0;
  return;
}
