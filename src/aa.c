/*
 * Anderson acceleration.
 *
 * x: input iterate
 * x_prev: previous input iterate
 * f: f(x) output of map f applied to x
 * g: x - f (error)
 * g_prev: previous error
 * s: x - x_prev
 * y: g - g_prev
 * d: s - y = f - f_prev
 *
 * capital letters are the variables stacked columnwise
 * idx tracks current index where latest quantities written
 * idx cycles from left to right columns in matrix
 *
 * Type-I:
 * return f = f - (S - Y) * ( S'Y + r I)^{-1} ( S'g )
 *
 * Type-II:
 * return f = f - (S - Y) * ( Y'Y + r I)^{-1} ( Y'g )
 *
 */

#include "aa.h"
#include "scs_blas.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define FILL_MEMORY_BEFORE_SOLVE (1)

#ifndef USE_LAPACK

typedef void *ACCEL_WORK;

AaWork *aa_init(aa_int dim, aa_int mem, aa_int type1, aa_float regularization,
                aa_float relaxation, aa_float safeguard_factor,
                aa_float max_weight_norm, aa_int verbosity) {
  return SCS_NULL;
}
aa_float aa_apply(aa_float *f, const aa_float *x, AaWork *a) {
  return 0;
}
aa_int aa_safeguard(aa_float *f_new, aa_float *x_new, AaWork *a) {
  return 0;
}
void aa_finish(AaWork *a) {
}
void aa_reset(AaWork *a) {
}

#else

#if PROFILING > 0

#define TIME_TIC                                                               \
  timer __t;                                                                   \
  tic(&__t);
#define TIME_TOC toc(__func__, &__t);

#include <time.h>
typedef struct timer {
  struct timespec tic;
  struct timespec toc;
} timer;

void tic(timer *t) {
  clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

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
  printf("%s - time: %8.4f milli-seconds.\n", str, time);
  return time;
}

#else

#define TIME_TIC
#define TIME_TOC

#endif

/* This file uses Anderson acceleration to improve the convergence of
 * a fixed point mapping.
 * At each iteration we need to solve a (small) linear system, we
 * do this using LAPACK ?gesv.
 */

/* contains the necessary parameters to perform aa at each step */
struct ACCEL_WORK {
  aa_int interval;
  aa_int dim;       /* variable dimension */
  aa_int iter;      /* current iteration */
  aa_int verbosity; /* verbosity level, 0 is no printing */
  aa_float *x;      /* average */
};

/*
 * API functions below this line, see aa.h for descriptions.
 */
AaWork *aa_init(aa_int dim, aa_int mem, aa_int type1, aa_float regularization,
                aa_float relaxation, aa_float safeguard_factor,
                aa_float max_weight_norm, aa_int verbosity) {
  TIME_TIC
  AaWork *a = (AaWork *)calloc(1, sizeof(AaWork));
  if (!a) {
    printf("Failed to allocate memory for AA.\n");
    return (AaWork *)0;
  }
  a->dim = dim;
  a->interval = mem; /* restart interval */
  a->verbosity = verbosity;
  a->x = (aa_float *)calloc(a->dim, sizeof(aa_float));
  aa_reset(a); /* initialize */
  TIME_TOC
  return a;
}

aa_float aa_apply(aa_float *f, const aa_float *x, AaWork *a) {
  TIME_TIC
  aa_int i;

  for (i = 0; i < a->dim; ++i) {
    /* accumulate the average */
    a->x[i] = (a->x[i] * a->iter + f[i]) / (a->iter + 1);
  }

  a->iter++;

  if (a->iter % a->interval == 0) {
    /* reset to average */
    if (a->verbosity > 0) {
      printf("reset to average iterate.\n");
    }
    for (i = 0; i < a->dim; ++i) {
      f[i] = a->x[i];
    }
    aa_reset(a);
  }

  TIME_TOC
  return 0.;
}

/* disabled */
aa_int aa_safeguard(aa_float *f_new, aa_float *x_new, AaWork *a) {
  return 0;
}

void aa_finish(AaWork *a) {
  if (a) {
    free(a->x);
    free(a);
  }
}

void aa_reset(AaWork *a) {
  /* to reset we simply set a->iter = 0 */
  if (a->verbosity > 0) {
    printf("AA reset.\n");
  }
  memset(a->x, 0, a->dim * sizeof(aa_float));
  a->iter = 0;
}

#endif
