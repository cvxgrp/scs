#include "accel.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h"
#include "util.h"

/* This file uses acceleration to improve the convergence of the ADMM iteration
 * z^+ = \phi(z). At each iteration we need to solve a (small) linear system, we
 * do this using LAPACK ?gesv.  If this fails then we just don't do any
 * acceleration this iteration, however we could fall back further to ?gelsy or
 * other more robust methods if we wanted to.
 */

#define MAX_ACCEL_PARAM_NORM (10.0)

struct SCS_ACCEL_WORK {
#ifdef USE_LAPACK
  scs_float *d_f;
  scs_float *d_g;
  scs_float *d_x;
  scs_float *delta_f;
  scs_float *delta_x;
  scs_float *f;
  scs_float *g;
  scs_float *x;
  scs_float *sol;
  scs_float *scratch;
  scs_float *mat;
  scs_float *tmp;
  blas_int *ipiv;
  scs_int k, l;
#endif
  scs_float total_accel_time;
};

#ifdef USE_LAPACK
void BLAS(gemv)(const char *trans, const blas_int *m, const blas_int *n,
                const scs_float *alpha, const scs_float *a, const blas_int *lda,
                const scs_float *x, const blas_int *incx, const scs_float *beta,
                scs_float *y, const blas_int *incy);
void BLAS(gesv)(blas_int *n, blas_int *nrhs, scs_float *a, blas_int *lda,
                blas_int *ipiv, scs_float *b, blas_int *ldb, blas_int *info);

static void update_mat(ScsAccelWork *a, scs_int idx) {
  /* use sol as scratch workspace here */
  DEBUG_FUNC
  scs_int i;
  scs_float *wrk = a->sol;
  scs_float *d_f = a->d_f;
  scs_float *d_x = a->d_x;
  scs_float *delta_x = a->delta_x;
  scs_float *delta_f = a->delta_f;
  scs_float *mat = a->mat;
  scs_int l = a->l;
  scs_int k = a->k;
  /* blas vars */
  blas_int twol = (blas_int)(2 * l);
  blas_int one = 1;
  blas_int bk = (blas_int)a->k;
  scs_float onef = 1.0;
  scs_float zerof = 0.0;

  scs_float ip = SCS(dot)(delta_x, delta_f, 2 * l);
  BLAS(gemv)
  ("Trans", &twol, &bk, &onef, d_x, &twol, delta_f, &one, &zerof, wrk, &one);
  SCS(add_scaled_array)(&(mat[idx * k]), wrk, k, -1.0);
  BLAS(gemv)
  ("Trans", &twol, &bk, &onef, d_f, &twol, delta_x, &one, &zerof, wrk, &one);
  for (i = 0; i < k; ++i) {
    mat[i * k + idx] -= wrk[i];
  }
  mat[idx * k + idx] += ip;
  RETURN;
}

static void update_accel_params(ScsWork *w, scs_int idx) {
  DEBUG_FUNC
  scs_float *d_f = w->accel->d_f;
  scs_float *d_g = w->accel->d_g;
  scs_float *d_x = w->accel->d_x;
  scs_float *f = w->accel->f;
  scs_float *g = w->accel->g;
  scs_float *x = w->accel->x;
  scs_float *delta_x = w->accel->delta_x;
  scs_float *delta_f = w->accel->delta_f;
  scs_int l = w->accel->l;

  /* copy g_prev into idx col of d_g */
  memcpy(&(d_g[idx * 2 * l]), g, sizeof(scs_float) * 2 * l);

  /* copy old col d_f into delta_f */
  memcpy(delta_f, &(d_f[idx * 2 * l]), sizeof(scs_float) * 2 * l);
  /* copy old col d_x into delta_x */
  memcpy(delta_x, &(d_x[idx * 2 * l]), sizeof(scs_float) * 2 * l);
  /* delta_f -= f_prev */
  SCS(add_scaled_array)(delta_f, f, 2 * l, -1.0);
  /* delta_x -= x_prev */
  SCS(add_scaled_array)(delta_x, x, 2 * l, -1.0);

  /* g = [u;v] */
  memcpy(g, w->u, sizeof(scs_float) * l);
  memcpy(&(g[l]), w->v, sizeof(scs_float) * l);
  /* x = [u_prev;v_prev] */
  memcpy(x, w->u_prev, sizeof(scs_float) * l);
  memcpy(&(x[l]), w->v_prev, sizeof(scs_float) * l);
  /* calculcate f = g - x */
  memcpy(f, g, sizeof(scs_float) * 2 * l);
  SCS(add_scaled_array)(f, x, 2 * l, -1.0);

  /* delta_f += f */
  SCS(add_scaled_array)(delta_f, f, 2 * l, 1.0);
  /* delta_x += x */
  SCS(add_scaled_array)(delta_x, x, 2 * l, 1.0);

  /* update mat = d_x'*d_f using delta_x, delta_f */
  update_mat(w->accel, idx);

  /* idx col of d_g = g_prev - g */
  SCS(add_scaled_array)(&(d_g[idx * 2 * l]), g, 2 * l, -1);
  /* idx col of d_f -= delta_f */
  SCS(add_scaled_array)(&(d_f[idx * 2 * l]), delta_f, 2 * l, -1);
  /* idx col of d_x -= delta_x */
  SCS(add_scaled_array)(&(d_x[idx * 2 * l]), delta_x, 2 * l, -1);
  RETURN;
}

ScsAccelWork *SCS(init_accel)(ScsWork *w) {
  DEBUG_FUNC
  ScsAccelWork *a = (ScsAccelWork *)scs_calloc(1, sizeof(ScsAccelWork));
  if (!a) {
    RETURN SCS_NULL;
  }
  a->l = w->m + w->n + 1;
  /* Use MIN to prevent not full rank matrices */
  a->k = MIN(w->n, w->stgs->acceleration_lookback);
  if (a->k <= 0) {
    RETURN a;
  }
  a->d_f = (scs_float *)scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->d_g = (scs_float *)scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->d_x = (scs_float *)scs_calloc(2 * a->l * a->k, sizeof(scs_float));
  a->f = (scs_float *)scs_calloc(2 * a->l, sizeof(scs_float));
  a->g = (scs_float *)scs_calloc(2 * a->l, sizeof(scs_float));
  a->x = (scs_float *)scs_calloc(2 * a->l, sizeof(scs_float));
  a->delta_f = (scs_float *)scs_calloc(2 * a->l, sizeof(scs_float));
  a->delta_x = (scs_float *)scs_calloc(2 * a->l, sizeof(scs_float));
  a->mat = (scs_float *)scs_calloc(a->k * a->k, sizeof(scs_float));
  a->tmp = (scs_float *)scs_calloc(a->k * a->k, sizeof(scs_float));
  a->sol = (scs_float *)scs_malloc(sizeof(scs_float) * 2 * a->l);
  a->scratch = (scs_float *)scs_malloc(sizeof(scs_float) * a->k);
  a->ipiv = (blas_int *)scs_malloc(sizeof(blas_int) * a->k);
  a->total_accel_time = 0.0;
  if (!a->d_f || !a->d_g || !a->f || !a->g || !a->scratch || !a->sol ||
      !a->d_x || !a->x || !a->scratch || !a->ipiv || !a->mat) {
    SCS(free_accel)(a);
    a = SCS_NULL;
  }
  RETURN a;
}

static scs_int solve_with_gesv(ScsAccelWork *a, scs_int len) {
  DEBUG_FUNC
  blas_int info = -1;
  blas_int twol = (blas_int)(2 * a->l);
  blas_int one = 1;
  blas_int blen = (blas_int)len;
  blas_int bk = (blas_int)a->k;
  scs_float neg_onef = -1.0;
  scs_float onef = 1.0;
  scs_float zerof = 0.0;
  scs_float *d_x = a->d_x;
  scs_float *mat = a->mat;
  scs_float *tmp = a->tmp;
  /* scratch = dX' f */
  BLAS(gemv)
  ("Trans", &twol, &blen, &onef, d_x, &twol, a->f, &one, &zerof, a->scratch,
   &one);
  /* copy mat into tmp since matrix is destroyed by gesv */
  memcpy(tmp, mat, a->k * a->k * sizeof(scs_float));
  /* scratch = (dX'dF) \ dX' f */
  BLAS(gesv)(&blen, &one, tmp, &bk, a->ipiv, a->scratch, &blen, &info);
  /* sol = g */
  memcpy(a->sol, a->g, sizeof(scs_float) * 2 * a->l);
  /* sol = sol - dG * scratch */
  BLAS(gemv)
  ("NoTrans", &twol, &blen, &neg_onef, a->d_g, &twol, a->scratch, &one, &onef,
   a->sol, &one);
  #if EXTRA_VERBOSE > 0
  scs_printf("norm of alphas %f\n", SCS(norm)(a->scratch, len));
  #endif
  if (SCS(norm)(a->scratch, len) >= MAX_ACCEL_PARAM_NORM) {
    RETURN -999;
  }
  RETURN(scs_int) info;
}

scs_int SCS(accelerate)(ScsWork *w, scs_int iter) {
  DEBUG_FUNC
  scs_int l = w->accel->l;
  scs_int k = w->accel->k;
  scs_int info = -1;
  SCS(timer) accel_timer;
  if (k <= 0) {
    RETURN 0;
  }
  SCS(tic)(&accel_timer);
  /* update df, d_g, d_x, f, g, x */
  update_accel_params(w, iter % k);
  if (iter == 0) {
    RETURN 0;
  }
  /* solve linear system, new point stored in sol */
  info = solve_with_gesv(w->accel, MIN(iter, k));
  w->accel->total_accel_time += SCS(tocq)(&accel_timer);
  /* check that info == 0 and fallback otherwise */
  if (info != 0) {
  #if EXTRA_VERBOSE > 0
    scs_printf(
        "Call to SCS(accelerate) failed with code %li, falling back to using "
        "no acceleration\n",
        (long)info);
  #endif
    RETURN 0;
  }
  /* set [u;v] = sol */
  memcpy(w->u, w->accel->sol, sizeof(scs_float) * l);
  memcpy(w->v, &(w->accel->sol[l]), sizeof(scs_float) * l);
  RETURN info;
}

void SCS(free_accel)(ScsAccelWork *a) {
  DEBUG_FUNC
  if (a) {
    if (a->d_f) {
      scs_free(a->d_f);
    }
    if (a->d_g) {
      scs_free(a->d_g);
    }
    if (a->d_x) {
      scs_free(a->d_x);
    }
    if (a->f) {
      scs_free(a->f);
    }
    if (a->g) {
      scs_free(a->g);
    }
    if (a->x) {
      scs_free(a->x);
    }
    if (a->delta_f) {
      scs_free(a->delta_f);
    }
    if (a->delta_x) {
      scs_free(a->delta_x);
    }
    if (a->sol) {
      scs_free(a->sol);
    }
    if (a->scratch) {
      scs_free(a->scratch);
    }
    if (a->mat) {
      scs_free(a->mat);
    }
    if (a->tmp) {
      scs_free(a->tmp);
    }
    if (a->ipiv) {
      scs_free(a->ipiv);
    }
    scs_free(a);
  }
  RETURN;
}

#else

ScsAccelWork *SCS(init_accel)(ScsWork *w) {
  ScsAccelWork *a = (ScsAccelWork *)scs_malloc(sizeof(ScsAccelWork));
  a->total_accel_time = 0.0;
  RETURN a;
}

void SCS(free_accel)(ScsAccelWork *a) {
  if (a) {
    scs_free(a);
  }
}

scs_int SCS(accelerate)(ScsWork *w, scs_int iter) { RETURN 0; }
#endif

char *SCS(get_accel_summary)(const ScsInfo *info, ScsAccelWork *a) {
  DEBUG_FUNC
  char *str = (char *)scs_malloc(sizeof(char) * 64);
  sprintf(str, "\tAcceleration: avg step time: %1.2es\n",
          a->total_accel_time / (info->iter + 1) / 1e3);
  a->total_accel_time = 0.0;
  RETURN str;
}
