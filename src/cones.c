#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h" /* contains BLAS(X) macros and type info */
#include "util.h"

/*
 * Cross-platform Complex Type Handling
 * MSVC uses struct layout; GCC/Clang uses C99 _Complex
 */
#if defined(_MSC_VER)
typedef struct {
  double real, imag;
} scs_blas_cdouble;
typedef struct {
  float real, imag;
} scs_blas_cfloat;
#ifndef SFLOAT
#define SCS_BLAS_COMPLEX_TYPE scs_blas_cdouble
#define SCS_BLAS_COMPLEX_CAST(x) ((scs_blas_cdouble *)(x))
#else
#define SCS_BLAS_COMPLEX_TYPE scs_blas_cfloat
#define SCS_BLAS_COMPLEX_CAST(x) ((scs_blas_cfloat *)(x))
#endif
#else
#include <complex.h>
#ifndef SFLOAT
#define SCS_BLAS_COMPLEX_CAST(x) ((double _Complex *)(x))
#define SCS_BLAS_COMPLEX_TYPE double _Complex
#else
#define SCS_BLAS_COMPLEX_CAST(x) ((float _Complex *)(x))
#define SCS_BLAS_COMPLEX_TYPE float _Complex
#endif
#endif

/* Constants */
#define BOX_CONE_MAX_ITERS (25)
#define POW_CONE_TOL (1e-9)
#define POW_CONE_MAX_ITERS (20)

/* Box cone limits (+ or -) taken to be INF */
#define MAX_BOX_VAL (1e15)

#ifdef USE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif

/* LAPACK / BLAS Function Prototypes */
void BLAS(syevr)(const char *jobz, const char *range, const char *uplo,
                 blas_int *n, scs_float *a, blas_int *lda, scs_float *vl,
                 scs_float *vu, blas_int *il, blas_int *iu, scs_float *abstol,
                 blas_int *m, scs_float *w, scs_float *z, blas_int *ldz,
                 blas_int *isuppz, scs_float *work, blas_int *lwork,
                 blas_int *iwork, blas_int *liwork, blas_int *info);

void BLASC(heevr)(const char *jobz, const char *range, const char *uplo,
                  blas_int *n, SCS_BLAS_COMPLEX_TYPE *a, blas_int *lda,
                  scs_float *vl, scs_float *vu, blas_int *il, blas_int *iu,
                  scs_float *abstol, blas_int *m, scs_float *w,
                  SCS_BLAS_COMPLEX_TYPE *z, blas_int *ldz, blas_int *isuppz,
                  SCS_BLAS_COMPLEX_TYPE *cwork, blas_int *lcwork,
                  scs_float *rwork, blas_int *lrwork, blas_int *iwork,
                  blas_int *liwork, blas_int *info);

blas_int BLAS(syrk)(const char *uplo, const char *trans, const blas_int *n,
                    const blas_int *k, const scs_float *alpha,
                    const scs_float *a, const blas_int *lda,
                    const scs_float *beta, scs_float *c, const blas_int *ldc);

blas_int BLASC(herk)(const char *uplo, const char *trans, const blas_int *n,
                     const blas_int *k, const scs_float *alpha,
                     const SCS_BLAS_COMPLEX_TYPE *a, const blas_int *lda,
                     const scs_float *beta, SCS_BLAS_COMPLEX_TYPE *c,
                     const blas_int *ldc);

void BLAS(scal)(const blas_int *n, const scs_float *sa, scs_float *sx,
                const blas_int *incx);
void BLASC(scal)(const blas_int *n, const SCS_BLAS_COMPLEX_TYPE *sa,
                 SCS_BLAS_COMPLEX_TYPE *sx, const blas_int *incx);

#ifdef USE_SPECTRAL_CONES
void BLAS(gesvd)(const char *jobu, const char *jobvt, const blas_int *m,
                 const blas_int *n, scs_float *a, const blas_int *lda,
                 scs_float *s, scs_float *u, const blas_int *ldu, scs_float *vt,
                 const blas_int *ldvt, scs_float *work, const blas_int *lwork,
                 blas_int *info);

/* Forward declare spectral matrix cone projections */
scs_int SCS(proj_logdet_cone)(scs_float *tvX, const scs_int n, ScsConeWork *c,
                              scs_int offset, bool *warmstart);
scs_int SCS(proj_nuclear_cone)(scs_float *tX, scs_int m, scs_int n,
                               ScsConeWork *c);
void SCS(proj_ell_one)(scs_float *tx, scs_int n, ScsConeWork *c);
scs_int SCS(proj_sum_largest_evals)(scs_float *tX, scs_int n, scs_int k,
                                    ScsConeWork *c);
#endif

#ifdef __cplusplus
}
#endif
#endif /* USE_LAPACK */

/* Forward declare exponential cone projection (exp_cone.c) */
scs_float SCS(proj_pd_exp_cone)(scs_float *v0, scs_int primal);

/*
 * Memory Management
 */

void SCS(free_cone)(ScsCone *k) {
  if (k) {
    if (k->bu)
      scs_free(k->bu);
    if (k->bl)
      scs_free(k->bl);
    if (k->q)
      scs_free(k->q);
    if (k->s)
      scs_free(k->s);
    if (k->cs)
      scs_free(k->cs);
    if (k->p)
      scs_free(k->p);
#ifdef USE_SPECTRAL_CONES
    if (k->d)
      scs_free(k->d);
    if (k->nuc_m)
      scs_free(k->nuc_m);
    if (k->nuc_n)
      scs_free(k->nuc_n);
    if (k->ell1)
      scs_free(k->ell1);
    if (k->sl_n)
      scs_free(k->sl_n);
    if (k->sl_k)
      scs_free(k->sl_k);
#endif
    scs_free(k);
  }
}

void SCS(deep_copy_cone)(ScsCone *dest, const ScsCone *src) {
  memcpy(dest, src, sizeof(ScsCone));

  /* Box cone */
  if (src->bsize > 1) {
    dest->bu = (scs_float *)scs_calloc(src->bsize - 1, sizeof(scs_float));
    dest->bl = (scs_float *)scs_calloc(src->bsize - 1, sizeof(scs_float));
    memcpy(dest->bu, src->bu, (src->bsize - 1) * sizeof(scs_float));
    memcpy(dest->bl, src->bl, (src->bsize - 1) * sizeof(scs_float));
  } else {
    dest->bu = SCS_NULL;
    dest->bl = SCS_NULL;
  }

  /* SOC */
  if (src->qsize > 0) {
    dest->q = (scs_int *)scs_calloc(src->qsize, sizeof(scs_int));
    memcpy(dest->q, src->q, src->qsize * sizeof(scs_int));
  } else {
    dest->q = SCS_NULL;
  }

  /* PSD */
  if (src->ssize > 0) {
    dest->s = (scs_int *)scs_calloc(src->ssize, sizeof(scs_int));
    memcpy(dest->s, src->s, src->ssize * sizeof(scs_int));
  } else {
    dest->s = SCS_NULL;
  }

  /* Complex PSD */
  if (src->cssize > 0) {
    dest->cs = (scs_int *)scs_calloc(src->cssize, sizeof(scs_int));
    memcpy(dest->cs, src->cs, src->cssize * sizeof(scs_int));
  } else {
    dest->cs = SCS_NULL;
  }

  /* Power */
  if (src->psize > 0) {
    dest->p = (scs_float *)scs_calloc(src->psize, sizeof(scs_float));
    memcpy(dest->p, src->p, src->psize * sizeof(scs_float));
  } else {
    dest->p = SCS_NULL;
  }

#ifdef USE_SPECTRAL_CONES
  /* Logdet */
  if (src->dsize > 0) {
    dest->d = (scs_int *)scs_calloc(src->dsize, sizeof(scs_int));
    memcpy(dest->d, src->d, src->dsize * sizeof(scs_int));
  } else {
    dest->d = SCS_NULL;
  }

  /* Nuclear */
  if (src->nucsize > 0) {
    dest->nuc_m = (scs_int *)scs_calloc(src->nucsize, sizeof(scs_int));
    dest->nuc_n = (scs_int *)scs_calloc(src->nucsize, sizeof(scs_int));
    memcpy(dest->nuc_m, src->nuc_m, src->nucsize * sizeof(scs_int));
    memcpy(dest->nuc_n, src->nuc_n, src->nucsize * sizeof(scs_int));
  } else {
    dest->nuc_m = SCS_NULL;
    dest->nuc_n = SCS_NULL;
  }

  /* Ell1 */
  if (src->ell1_size > 0) {
    dest->ell1 = (scs_int *)scs_calloc(src->ell1_size, sizeof(scs_int));
    memcpy(dest->ell1, src->ell1, src->ell1_size * sizeof(scs_int));
  } else {
    dest->ell1 = SCS_NULL;
  }

  /* Sum Largest Eigenvalues */
  if (src->sl_size > 0) {
    dest->sl_n = (scs_int *)scs_calloc(src->sl_size, sizeof(scs_int));
    dest->sl_k = (scs_int *)scs_calloc(src->sl_size, sizeof(scs_int));
    memcpy(dest->sl_n, src->sl_n, src->sl_size * sizeof(scs_int));
    memcpy(dest->sl_k, src->sl_k, src->sl_size * sizeof(scs_int));
  } else {
    dest->sl_n = SCS_NULL;
    dest->sl_k = SCS_NULL;
  }
#endif
}

/*
 * Helper Functions
 */

static inline scs_int get_sd_cone_size(scs_int s) {
  return (s * (s + 1)) / 2;
}
static inline scs_int get_csd_cone_size(scs_int cs) {
  return cs * cs;
}

void SCS(set_r_y)(const ScsConeWork *c, scs_float scale, scs_float *r_y) {
  scs_int i;
  /* z cone */
  for (i = 0; i < c->k->z; ++i) {
    /* set rho_y small for z, similar to rho_x term, since z corresponds to
     * dual free cone, this effectively decreases penalty on those entries
     * and lets them be determined almost entirely by the linear system solve
     */
    r_y[i] = 1.0 / (1000. * scale);
  }
  /* Remaining cones */
  for (i = c->k->z; i < c->m; ++i) {
    r_y[i] = 1.0 / scale;
  }
}

/* The function f aggregates the entries within each cone */
void SCS(enforce_cone_boundaries)(const ScsConeWork *c, scs_float *vec,
                                  scs_float (*f)(const scs_float *, scs_int)) {
  scs_int i, j, delta;
  scs_int count = c->cone_boundaries[0];
  scs_float wrk;
  for (i = 1; i < c->cone_boundaries_len; ++i) {
    delta = c->cone_boundaries[i];
    wrk = f(&(vec[count]), delta);
    for (j = count; j < count + delta; ++j) {
      vec[j] = wrk;
    }
    count += delta;
  }
}

/*
 * Boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size strictly
 * larger than 1, boundaries malloc-ed here so should be freed.
 */
void set_cone_boundaries(const ScsCone *k, ScsConeWork *c) {
  scs_int i, count = 0;
#ifdef USE_SPECTRAL_CONES
  scs_int total_cones = k->qsize + k->ssize + k->cssize + k->ed + k->ep +
                        k->psize + k->dsize + k->nucsize + k->ell1_size +
                        k->sl_size;
#else
  scs_int total_cones =
      k->qsize + k->ssize + k->cssize + k->ed + k->ep + k->psize;
#endif
  scs_int *b = (scs_int *)scs_calloc(total_cones + 1, sizeof(scs_int));

  /* Cones that can be scaled independently */
  b[count++] = k->z + k->l + k->bsize;
  for (i = 0; i < k->qsize; ++i)
    b[count++] = k->q[i];
  for (i = 0; i < k->ssize; ++i)
    b[count++] = get_sd_cone_size(k->s[i]);
  for (i = 0; i < k->cssize; ++i)
    b[count++] = get_csd_cone_size(k->cs[i]);
  for (i = 0; i < k->ep + k->ed; ++i)
    b[count++] = 3;
  for (i = 0; i < k->psize; ++i)
    b[count++] = 3;

#ifdef USE_SPECTRAL_CONES
  for (i = 0; i < k->dsize; ++i)
    b[count++] = get_sd_cone_size(k->d[i]) + 2;
  for (i = 0; i < k->nucsize; ++i)
    b[count++] = k->nuc_m[i] * k->nuc_n[i] + 1;
  for (i = 0; i < k->ell1_size; ++i)
    b[count++] = k->ell1[i] + 1;
  for (i = 0; i < k->sl_size; ++i)
    b[count++] = get_sd_cone_size(k->sl_n[i]) + 1;
#endif

  c->cone_boundaries = b;
  c->cone_boundaries_len = total_cones + 1;
}

static scs_int get_full_cone_dims(const ScsCone *k) {
  scs_int i, dims = k->z + k->l + k->bsize;
  for (i = 0; i < k->qsize; ++i)
    dims += k->q[i];
  for (i = 0; i < k->ssize; ++i)
    dims += get_sd_cone_size(k->s[i]);
  for (i = 0; i < k->cssize; ++i)
    dims += get_csd_cone_size(k->cs[i]);
  dims += 3 * (k->ed + k->ep + k->psize);
#ifdef USE_SPECTRAL_CONES
  for (i = 0; i < k->dsize; ++i)
    dims += get_sd_cone_size(k->d[i]) + 2;
  for (i = 0; i < k->nucsize; ++i)
    dims += k->nuc_m[i] * k->nuc_n[i] + 1;
  for (i = 0; i < k->ell1_size; ++i)
    dims += k->ell1[i] + 1;
  for (i = 0; i < k->sl_size; ++i)
    dims += get_sd_cone_size(k->sl_n[i]) + 1;
#endif
  return dims;
}

scs_int SCS(validate_cones)(const ScsData *d, const ScsCone *k) {
  scs_int i;
  if (get_full_cone_dims(k) != d->m) {
    scs_printf("Error: Cone dims %li != rows in A %li\n",
               (long)get_full_cone_dims(k), (long)d->m);
    return -1;
  }
  if (k->z && k->z < 0) {
    scs_printf("free cone dimension error\n");
    return -1;
  }
  if (k->l && k->l < 0) {
    scs_printf("lp cone dimension error\n");
    return -1;
  }
  if (k->bsize) {
    if (k->bsize < 0) {
      scs_printf("box cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->bsize - 1; ++i) {
      if (k->bl[i] > k->bu[i]) {
        scs_printf("infeasible: box lower bound larger than upper bound\n");
        return -1;
      }
    }
  }
  if (k->qsize && k->q) {
    if (k->qsize < 0) {
      scs_printf("soc cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->qsize; ++i) {
      if (k->q[i] < 0) {
        scs_printf("soc cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->ssize && k->s) {
    if (k->ssize < 0) {
      scs_printf("sd cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->ssize; ++i) {
      if (k->s[i] < 0) {
        scs_printf("sd cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->cssize && k->cs) {
    if (k->cssize < 0) {
      scs_printf("complex psd cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->cssize; ++i) {
      if (k->cs[i] < 0) {
        scs_printf("complex psd cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->ed && k->ed < 0) {
    scs_printf("ep cone dimension error\n");
    return -1;
  }
  if (k->ep && k->ep < 0) {
    scs_printf("ed cone dimension error\n");
    return -1;
  }
  if (k->psize && k->p) {
    if (k->psize < 0) {
      scs_printf("power cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->psize; ++i) {
      if (k->p[i] < -1 || k->p[i] > 1) {
        scs_printf("power cone error, values must be in [-1,1]\n");
        return -1;
      }
    }
  }
#ifdef USE_SPECTRAL_CONES
  if (k->dsize && k->d) {
    if (k->dsize < 0) {
      scs_printf("logdet cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->dsize; ++i) {
      if (k->d[i] < 1) {
        scs_printf("logdet cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->nucsize && k->nuc_m && k->nuc_n) {
    if (k->nucsize < 0) {
      scs_printf("nuclear cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->nucsize; ++i) {
      if (k->nuc_m[i] < 1 || k->nuc_n[i] < 1 || k->nuc_n[i] > k->nuc_m[i]) {
        scs_printf("nuclear norm cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->ell1_size && k->ell1) {
    if (k->ell1_size < 0) {
      scs_printf("ell1 cone dimension error\n");
      return -1;
    }
    for (i = 0; i < k->ell1_size; ++i) {
      if (k->ell1[i] < 1) {
        scs_printf("ell1 cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->sl_size && k->sl_n && k->sl_k) {
    for (i = 0; i < k->sl_size; ++i) {
      if ((k->sl_k[i] >= k->sl_n[i]) || k->sl_k[i] <= 0) {
        scs_printf("sum-of-largest-eigenvalues cone dimension error\n");
        return -1;
      }
    }
  }
#endif
  return 0;
}

void SCS(finish_cone)(ScsConeWork *c) {
  if (!c)
    return;
#ifdef USE_LAPACK
  if (c->Xs)
    scs_free(c->Xs);
  if (c->cXs)
    scs_free(c->cXs);
  if (c->Z)
    scs_free(c->Z);
  if (c->cZ)
    scs_free(c->cZ);
  if (c->e)
    scs_free(c->e);
  if (c->isuppz)
    scs_free(c->isuppz);
  if (c->work)
    scs_free(c->work);
  if (c->iwork)
    scs_free(c->iwork);
  if (c->cwork)
    scs_free(c->cwork);
  /* c->rwork is aliased to c->work in setup, no free needed */
#endif
  if (c->cone_boundaries)
    scs_free(c->cone_boundaries);
  if (c->s)
    scs_free(c->s);

#ifdef USE_SPECTRAL_CONES
  if (c->work_logdet)
    scs_free(c->work_logdet);
  if (c->saved_log_projs)
    scs_free(c->saved_log_projs);
  if (c->s_nuc)
    scs_free(c->s_nuc);
  if (c->u_nuc)
    scs_free(c->u_nuc);
  if (c->vt_nuc)
    scs_free(c->vt_nuc);
  if (c->work_nuc)
    scs_free(c->work_nuc);
  if (c->work_sum_of_largest)
    scs_free(c->work_sum_of_largest);
  if (c->log_cone_warmstarts)
    scs_free(c->log_cone_warmstarts);
  if (c->work_ell1)
    scs_free(c->work_ell1);
  if (c->work_ell1_proj)
    scs_free(c->work_ell1_proj);
#endif
  scs_free(c);
}

char *SCS(get_cone_header)(const ScsCone *k) {
  char *tmp = (char *)scs_malloc(512);
  scs_int i, count;

  sprintf(tmp, "cones: ");
  if (k->z)
    sprintf(tmp + strlen(tmp), "\t  z: primal zero / dual free vars: %li\n",
            (long)k->z);
  if (k->l)
    sprintf(tmp + strlen(tmp), "\t  l: linear vars: %li\n", (long)k->l);
  if (k->bsize)
    sprintf(tmp + strlen(tmp), "\t  b: box cone vars: %li\n", (long)k->bsize);

  if (k->qsize) {
    count = 0;
    for (i = 0; i < k->qsize; ++i)
      count += k->q[i];
    sprintf(tmp + strlen(tmp), "\t  q: soc vars: %li, qsize: %li\n",
            (long)count, (long)k->qsize);
  }
  if (k->ssize) {
    count = 0;
    for (i = 0; i < k->ssize; ++i)
      count += get_sd_cone_size(k->s[i]);
    sprintf(tmp + strlen(tmp), "\t  s: psd vars: %li, ssize: %li\n",
            (long)count, (long)k->ssize);
  }
  if (k->cssize) {
    count = 0;
    for (i = 0; i < k->cssize; ++i)
      count += get_csd_cone_size(k->cs[i]);
    sprintf(tmp + strlen(tmp), "\t  cs: complex psd vars: %li, cssize: %li\n",
            (long)count, (long)k->cssize);
  }
  if (k->ep || k->ed) {
    sprintf(tmp + strlen(tmp), "\t  e: exp vars: %li, dual exp vars: %li\n",
            (long)(3 * k->ep), (long)(3 * k->ed));
  }
  if (k->psize) {
    sprintf(tmp + strlen(tmp), "\t  p: primal + dual power vars: %li\n",
            (long)(3 * k->psize));
  }
#ifdef USE_SPECTRAL_CONES
  scs_int ell1_vars, log_vars, nuc_vars, sl_vars;
  log_vars = 0;
  if (k->dsize && k->d) {
    for (i = 0; i < k->dsize; i++) {
      log_vars += get_sd_cone_size(k->d[i]) + 2;
    }
    sprintf(tmp + strlen(tmp), "\t  d: logdet vars: %li, dsize: %li\n",
            (long)log_vars, (long)k->dsize);
  }
  nuc_vars = 0;
  if (k->nucsize && k->nuc_m && k->nuc_n) {
    for (i = 0; i < k->nucsize; i++) {
      nuc_vars += k->nuc_m[i] * k->nuc_n[i] + 1;
    }
    sprintf(tmp + strlen(tmp), "\t  nuc: nuclear vars: %li, nucsize: %li\n",
            (long)nuc_vars, (long)k->nucsize);
  }
  ell1_vars = 0;
  if (k->ell1_size && k->ell1) {
    for (i = 0; i < k->ell1_size; ++i) {
      ell1_vars += k->ell1[i];
    }
    sprintf(tmp + strlen(tmp), "\t  ell1: ell1 vars: %li, ell1_size: %li\n",
            (long)ell1_vars, (long)k->ell1_size);
  }

  sl_vars = 0;
  if (k->sl_size && k->sl_n) {
    for (i = 0; i < k->sl_size; ++i) {
      sl_vars += get_sd_cone_size(k->sl_n[i]) + 1;
    }
    sprintf(tmp + strlen(tmp), "\t  sl: sl vars: %li, sl_size: %li\n",
            (long)sl_vars, (long)k->sl_size);
  }
#endif
  return tmp;
}

/*
 * Workspace Setup
 * Consolidated setup for Real PSD, Complex PSD, and Spectral cones.
 */
static scs_int set_up_cone_work_spaces(ScsConeWork *c, const ScsCone *k) {
  scs_int i;
#ifdef USE_LAPACK
  /* Max dim for eigenvalues (e) and integer work (isuppz) */
  blas_int n_max = 1;
  blas_int n_max_real = 1; /* Max dim for Real PSD matrix (Xs) */
  blas_int n_max_csd = 1;  /* Max dim for Complex PSD matrix (cXs) */

  /* LAPACK Query Variables */
  blas_int neg_one = -1, info = 0, m = 0;
  blas_int d_i = 0;
  scs_float d_f = 0.0, abstol = -1.0;
  scs_float wkopt = 0.0;
  blas_int iwkopt = 0;
  scs_complex_float lcwork_opt = {0.0};
  scs_float lrwork_opt = 0.0;
  blas_int liwork_opt_c = 0;

  /* Max workspace sizes across all cone types */
  blas_int lwork_max = 0;
  blas_int liwork_max = 0;

  /* Calculate Dimensions */
  for (i = 0; i < k->ssize; ++i) {
    n_max = MAX(n_max, (blas_int)k->s[i]);
    n_max_real = MAX(n_max_real, (blas_int)k->s[i]);
  }
  for (i = 0; i < k->cssize; ++i) {
    n_max = MAX(n_max, (blas_int)k->cs[i]);
    n_max_csd = MAX(n_max_csd, (blas_int)k->cs[i]);
  }

#ifdef USE_SPECTRAL_CONES
  blas_int n_max_logdet = 1;
  blas_int n_logdet_total = 0;
  blas_int n_max_sl = 1;

  c->log_cone_warmstarts = (bool *)scs_calloc(k->dsize, sizeof(bool));

  /* Spectral Dimensions */
  for (i = 0; i < k->sl_size; ++i)
    n_max_sl = MAX(n_max_sl, (blas_int)k->sl_n[i]);
  for (i = 0; i < k->dsize; ++i) {
    n_logdet_total += (blas_int)k->d[i];
    n_max_logdet = MAX(n_max_logdet, (blas_int)k->d[i]);
  }

  /* Logdet Allocation */
  if (k->dsize > 0) {
    c->work_logdet =
        (scs_float *)scs_calloc(22 * n_max_logdet + 122, sizeof(scs_float));
    c->saved_log_projs = (scs_float *)scs_calloc(2 * k->dsize + n_logdet_total,
                                                 sizeof(scs_float));
    if (!c->work_logdet || !c->saved_log_projs)
      return -1;
  }
  /* Sum Largest Allocation */
  if (k->sl_size > 0) {
    c->work_sum_of_largest =
        (scs_float *)scs_calloc(n_max_sl * n_max_sl, sizeof(scs_float));
    if (!c->work_sum_of_largest)
      return -1;
  }

  /* Update Max Real Dims based on Spectral needs */
  n_max = MAX(n_max, n_max_logdet);
  n_max = MAX(n_max, n_max_sl);
  n_max_real = MAX(n_max_real, n_max_logdet);
  n_max_real = MAX(n_max_real, n_max_sl);
#endif

  /* Allocate standard eigenvalue buffers
   * 'e' stores eigenvalues, 'isuppz' supports them. Shared by Real/Complex.
   */
  c->e = (scs_float *)scs_calloc(n_max, sizeof(scs_float));
  c->isuppz = (blas_int *)scs_calloc(MAX(2, 2 * n_max), sizeof(blas_int));
  if (!c->e || !c->isuppz)
    return -1;

  /* 1. Real PSD Workspace Query (syevr) */
  if (k->ssize > 0
#ifdef USE_SPECTRAL_CONES
      || k->dsize > 0 || k->sl_size > 0
#endif
  ) {
    c->Xs = (scs_float *)scs_calloc(n_max_real * n_max_real, sizeof(scs_float));
    c->Z = (scs_float *)scs_calloc(n_max_real * n_max_real, sizeof(scs_float));
    if (!c->Xs || !c->Z)
      return -1;

    BLAS(syevr)("V", "A", "L", &n_max_real, c->Xs, &n_max_real, &d_f, &d_f,
                &d_i, &d_i, &abstol, &m, c->e, c->Z, &n_max_real, c->isuppz,
                &wkopt, &neg_one, &iwkopt, &neg_one, &info);

    if (info != 0) {
      scs_printf("FATAL: syevr workspace query failure, info = %li\n",
                 (long)info);
      return -1;
    }
    lwork_max = MAX(lwork_max, (blas_int)(wkopt + 1));
    liwork_max = MAX(liwork_max, iwkopt);
  }

  /* 2. Complex PSD Workspace Query (heevr) */
  if (k->cssize > 0) {
    c->cXs = (scs_complex_float *)scs_calloc(n_max_csd * n_max_csd,
                                             sizeof(scs_complex_float));
    c->cZ = (scs_complex_float *)scs_calloc(n_max_csd * n_max_csd,
                                            sizeof(scs_complex_float));
    if (!c->cXs || !c->cZ)
      return -1;

    BLASC(heevr)("V", "A", "L", &n_max_csd, SCS_BLAS_COMPLEX_CAST(c->cXs),
                 &n_max_csd, &d_f, &d_f, &d_i, &d_i, &abstol, &m, c->e,
                 SCS_BLAS_COMPLEX_CAST(c->cZ), &n_max_csd, c->isuppz,
                 SCS_BLAS_COMPLEX_CAST(&lcwork_opt), &neg_one, &lrwork_opt,
                 &neg_one, &liwork_opt_c, &neg_one, &info);

    if (info != 0) {
      scs_printf("FATAL: heevr workspace query failure, info = %li\n",
                 (long)info);
      return -1;
    }

    c->lcwork = (blas_int)(lcwork_opt[0]);
    c->cwork =
        (scs_complex_float *)scs_calloc(c->lcwork, sizeof(scs_complex_float));
    if (!c->cwork)
      return -1;

    /* heevr uses a real 'rwork' array. We alias this to the shared real 'work'
     * array */
    lwork_max = MAX(lwork_max, (blas_int)lrwork_opt);
    liwork_max = MAX(liwork_max, liwork_opt_c);
  }

#ifdef USE_SPECTRAL_CONES
  /* 3. Nuclear Norm Workspace (gesvd) */
  if (k->nucsize > 0) {
    blas_int m_max = 1, n_max_n = 1;
    for (i = 0; i < k->nucsize; ++i) {
      m_max = MAX(m_max, (blas_int)k->nuc_m[i]);
      n_max_n = MAX(n_max_n, (blas_int)k->nuc_n[i]);
    }

    c->s_nuc = (scs_float *)scs_calloc(n_max_n, sizeof(scs_float));
    c->u_nuc = (scs_float *)scs_calloc(m_max * n_max_n, sizeof(scs_float));
    c->vt_nuc = (scs_float *)scs_calloc(n_max_n * n_max_n, sizeof(scs_float));
    if (!c->s_nuc || !c->u_nuc || !c->vt_nuc)
      return -1;

    BLAS(gesvd)("S", "A", &m_max, &n_max_n, c->u_nuc, &m_max, c->s_nuc,
                c->u_nuc, &m_max, c->vt_nuc, &n_max_n, &wkopt, &neg_one, &info);

    if (info != 0)
      return -1;

    /* We allocate work_nuc separately to avoid aliasing risks during SVD
     * perations */
    c->lwork_nuc = (blas_int)(wkopt + 1);
    c->work_nuc = (scs_float *)scs_calloc(c->lwork_nuc, sizeof(scs_float));
    if (!c->work_nuc)
      return -1;
  }
#endif

  /* Final Consolidated Allocation
   * c->work aliases 'work' for syevr and 'rwork' for heevr
   */
  if (lwork_max > 0) {
    c->lwork = lwork_max;
    c->work = (scs_float *)scs_calloc(c->lwork, sizeof(scs_float));
    if (!c->work)
      return -1;
  }
  if (liwork_max > 0) {
    c->liwork = liwork_max;
    c->iwork = (blas_int *)scs_calloc(c->liwork, sizeof(blas_int));
    if (!c->iwork)
      return -1;
  }

  return 0;
#else
  /* Non-LAPACK fallback check */
  if (k->ssize > 0 || k->cssize > 0) {
    for (i = 0; i < k->ssize; i++) {
      if (k->s[i] > 1) {
        scs_printf("FATAL: SDP/Complex SDP requires BLAS/LAPACK.\n");
        return -1;
      }
    }
    for (i = 0; i < k->cssize; i++) {
      if (k->cs[i] > 1) {
        scs_printf("FATAL: SDP/Complex SDP requires BLAS/LAPACK.\n");
        return -1;
      }
    }
  }
#ifdef USE_SPECTRAL_CONES
  if (k->dsize > 0 || k->nucsize > 0 || k->sl_size > 0) {
    scs_printf("FATAL: Spectral cones require BLAS/LAPACK.\n");
    return -1;
  }
#endif
  return 0;
#endif
}

#ifdef USE_SPECTRAL_CONES
static scs_int set_up_ell1_cone_work_space(ScsConeWork *c, const ScsCone *k) {
  scs_int i, n_max = 0;
  if (k->ell1_size > 0) {
    for (i = 0; i < k->ell1_size; ++i)
      n_max = MAX(k->ell1[i], n_max);
    c->work_ell1 = (Value_index *)scs_calloc(n_max, sizeof(Value_index));
    c->work_ell1_proj = (scs_float *)scs_calloc(n_max + 1, sizeof(scs_float));
    if (!c->work_ell1 || !c->work_ell1_proj)
      return -1;
  }
  return 0;
}
#endif

/*
 * Projection: Real Semi-Definite Cone
 */
static scs_int proj_semi_definite_cone(scs_float *X, const scs_int n,
                                       ScsConeWork *c) {
  if (n == 0)
    return 0;
  if (n == 1) {
    X[0] = MAX(X[0], 0.);
    return 0;
  }

#ifdef USE_LAPACK
  scs_int i;
  blas_int nb = (blas_int)n;
  blas_int nb_p1 = (blas_int)(n + 1);
  blas_int info = 0, one_int = 1, ncols_z = 0;
  scs_float zero = 0., one = 1., sqrt2 = SQRTF(2.0), sqrt2_inv = 1.0 / sqrt2;
  scs_float abstol = -1.0, d_f = 0.0, sq_eig;
  blas_int m = 0, d_i = 0;
  scs_int first_idx = -1;

  /* Copy lower triangular part to full matrix buffer Xs */
  for (i = 0; i < n; ++i) {
    memcpy(&(c->Xs[i * (n + 1)]), &(X[i * n - ((i - 1) * i) / 2]),
           (n - i) * sizeof(scs_float));
  }

  /* Scale diagonals by sqrt(2) */
  BLAS(scal)(&nb, &sqrt2, c->Xs, &nb_p1);

  /* Solve Eigenproblem: Xs = Z * diag(e) * Z' */
  BLAS(syevr)("V", "A", "L", &nb, c->Xs, &nb, &d_f, &d_f, &d_i, &d_i, &abstol,
              &m, c->e, c->Z, &nb, c->isuppz, c->work, &c->lwork, c->iwork,
              &c->liwork, &info);
  if (info != 0)
    return (int)info;

  /* Filter negative eigenvalues and scale eigenvectors */
  /* Note: e is in ascending order */
  for (i = 0; i < n; ++i) {
    if (c->e[i] > 0) {
      if (first_idx == -1) {
        first_idx = i;
      }
      sq_eig = SQRTF(c->e[i]);
      BLAS(scal)(&nb, &sq_eig, &c->Z[i * n], &one_int);
    }
  }
  if (first_idx == -1) {
    /* there are no positive eigenvalues, set X to 0 and return */
    memset(X, 0, sizeof(scs_float) * get_sd_cone_size(n));
    return 0;
  }
  /* Reconstruct Xs = Z * Z' */
  ncols_z = (blas_int)(n - first_idx);
  BLAS(syrk)("Lower", "NoTrans", &nb, &ncols_z, &one, &c->Z[first_idx * n], &nb,
             &zero, c->Xs, &nb);

  /* Rescale diagonals by 1/sqrt(2) */
  BLAS(scal)(&nb, &sqrt2_inv, c->Xs, &nb_p1);

  /* Extract lower triangular matrix back to X */
  for (i = 0; i < n; ++i) {
    memcpy(&(X[i * n - ((i - 1) * i) / 2]), &(c->Xs[i * (n + 1)]),
           (n - i) * sizeof(scs_float));
  }
  return 0;
#else
  return -1;
#endif
}

/*
 * Projection: Complex Semi-Definite Cone
 */
static scs_int proj_complex_semi_definite_cone(scs_float *X, const scs_int n,
                                               ScsConeWork *c) {
  if (n == 0)
    return 0;
  if (n == 1) {
    X[0] = MAX(X[0], 0.);
    return 0;
  }

#ifdef USE_LAPACK
  scs_int i;
  blas_int nb = (blas_int)n;
  blas_int nb_p1 = (blas_int)(n + 1);
  blas_int info = 0, one_int = 1, d_i = 0, ncols_z = 0;
  scs_float zero = 0., one = 1., abstol = -1.0, d_f = 0.0;
  blas_int m = 0;
  scs_int first_idx = -1;

  /* Complex constants */
  scs_complex_float csqrt2 = {0.0}, csqrt2_inv = {0.0}, csq_eig = {0.0};
  csqrt2[0] = SQRTF(2.0);
  csqrt2_inv[0] = 1.0 / csqrt2[0];

  /* Unpack X (real array) into cXs (complex matrix) */
  for (i = 0; i < n - 1; ++i) {
    c->cXs[i * (n + 1)][0] = X[i * (2 * n - i)]; /* Diagonal (Real) */
    c->cXs[i * (n + 1)][1] = 0.0;
    memcpy(&(c->cXs[i * (n + 1) + 1]), &(X[i * (2 * n - i) + 1]),
           2 * (n - i - 1) * sizeof(scs_float));
  }
  c->cXs[n * n - 1][0] = X[n * n - 1]; /* Last element */
  c->cXs[n * n - 1][1] = 0.0;

  /* Scale diagonals by sqrt*/
  BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2),
              SCS_BLAS_COMPLEX_CAST(c->cXs), &nb_p1);

  /* Solve Eigenproblem. Note: c->work acts as rwork here. */
  BLASC(heevr)("V", "A", "L", &nb, SCS_BLAS_COMPLEX_CAST(c->cXs), &nb, &d_f,
               &d_f, &d_i, &d_i, &abstol, &m, c->e,
               SCS_BLAS_COMPLEX_CAST(c->cZ), &nb, c->isuppz,
               SCS_BLAS_COMPLEX_CAST(c->cwork), &c->lcwork, c->work, &c->lwork,
               c->iwork, &c->liwork, &info);
  if (info != 0)
    return (int)info;

  /* Reconstruct */
  for (i = 0; i < n; ++i) {
    if (c->e[i] > 0) {
      if (first_idx == -1) {
        first_idx = i;
      }
      csq_eig[0] = SQRTF(c->e[i]);
      BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csq_eig),
                  SCS_BLAS_COMPLEX_CAST(&c->cZ[i * n]), &one_int);
    }
  }
  if (first_idx == -1) {
    /* there are no positive eigenvalues, set X to 0 and return */
    memset(X, 0, sizeof(scs_float) * get_csd_cone_size(n));
    return 0;
  }

  /* cXs = cZ * cZ' */
  ncols_z = (blas_int)(n - first_idx);
  BLASC(herk)("Lower", "NoTrans", &nb, &ncols_z, &one,
              SCS_BLAS_COMPLEX_CAST(&c->cZ[first_idx * n]), &nb, &zero,
              SCS_BLAS_COMPLEX_CAST(c->cXs), &nb);

  /* Rescale diagonals */
  BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2_inv),
              SCS_BLAS_COMPLEX_CAST(c->cXs), &nb_p1);

  /* Repack into X */
  for (i = 0; i < n - 1; ++i) {
    X[i * (2 * n - i)] = c->cXs[i * (n + 1)][0];
    memcpy(&(X[i * (2 * n - i) + 1]), &(c->cXs[i * (n + 1) + 1]),
           2 * (n - i - 1) * sizeof(scs_float));
  }
  X[n * n - 1] = c->cXs[n * n - 1][0];
  return 0;
#else
  return -1;
#endif
}

/*
 * Projection: Box Cone
 */
static void normalize_box_cone(ScsCone *k, scs_float *D, scs_int bsize) {
  scs_int j;
  scs_float factor;
  for (j = 0; j < bsize - 1; j++) {
    factor = D ? D[j + 1] / D[0] : 1.0;

    if (k->bu[j] >= MAX_BOX_VAL)
      k->bu[j] = INFINITY;
    else
      k->bu[j] *= factor;

    if (k->bl[j] <= -MAX_BOX_VAL)
      k->bl[j] = -INFINITY;
    else
      k->bl[j] *= factor;
  }
}

/* Project onto { (t, s) | t * l <= s <= t * u, t >= 0 }, Newton's method on t
   tx = [t; s], total length = bsize, under Euclidean metric 1/r_box.
*/
static scs_float proj_box_cone(scs_float *tx, const scs_float *bl,
                               const scs_float *bu, scs_int bsize,
                               scs_float t_wm, scs_float *r_box) {
  scs_float *x = &(tx[1]);
  scs_float gt, ht, t = t_wm, t_prev, r;
  scs_float rho_t = 1.0;
  scs_float *rho = SCS_NULL;
  scs_int iter, j;

  if (bsize == 1) {
    tx[0] = MAX(tx[0], 0.0);
    return tx[0];
  }
  if (r_box) {
    rho_t = 1.0 / r_box[0];
    rho = &(r_box[1]);
  }

  /* Newton's method for t */
  for (iter = 0; iter < BOX_CONE_MAX_ITERS; iter++) {
    t_prev = t;
    gt = rho_t * (t - tx[0]);
    ht = rho_t;

    for (j = 0; j < bsize - 1; j++) {
      r = rho ? 1.0 / rho[j] : 1.0;
      if (x[j] > t * bu[j]) {
        gt += r * (t * bu[j] - x[j]) * bu[j];
        ht += r * bu[j] * bu[j];
      } else if (x[j] < t * bl[j]) {
        gt += r * (t * bl[j] - x[j]) * bl[j];
        ht += r * bl[j] * bl[j];
      }
    }

    t = MAX(t - gt / MAX(ht, 1e-8), 0.0);

    /*
     * TODO: sometimes this check can fail (ie, declare convergence before it
     * should) if ht is very large, which can happen with some pathological
     * problems.
     */
    if (ABS(gt / MAX(ht, 1e-6)) < 1e-12 * MAX(t, 1.) ||
        ABS(t - t_prev) < 1e-11 * MAX(t, 1.))
      break;
  }

  /* Apply calculated t to x */
  for (j = 0; j < bsize - 1; j++) {
    if (x[j] > t * bu[j])
      x[j] = t * bu[j];
    else if (x[j] < t * bl[j])
      x[j] = t * bl[j];
  }
  tx[0] = t;
  return t;
}

/*
 * Projection: Second Order Cone
 */
static void proj_soc(scs_float *x, scs_int q) {
  if (q <= 0)
    return;
  if (q == 1) {
    x[0] = MAX(x[0], 0.);
    return;
  }

  scs_float v1 = x[0];
  scs_float s = SCS(norm_2)(&(x[1]), q - 1);
  scs_float alpha = (s + v1) / 2.0;

  if (s <= v1)
    return;       /* Inside cone */
  if (s <= -v1) { /* Below dual cone */
    memset(x, 0, q * sizeof(scs_float));
  } else { /* Projection */
    x[0] = alpha;
    SCS(scale_array)(&(x[1]), alpha / s, q - 1);
  }
}

/*
 * Projection: Power Cone
 */
static scs_float pow_calc_x(scs_float r, scs_float xh, scs_float rh,
                            scs_float a) {
  scs_float x = 0.5 * (xh + SQRTF(xh * xh + 4 * a * (rh - r) * r));
  return MAX(x, 1e-12);
}

static scs_float pow_calc_fp(scs_float x, scs_float y, scs_float dxdr,
                             scs_float dydr, scs_float a) {
  return POWF(x, a) * POWF(y, (1 - a)) * (a * dxdr / x + (1 - a) * dydr / y) -
         1;
}

static void proj_power_cone(scs_float *v, scs_float a) {
  scs_float xh = v[0], yh = v[1], rh = ABS(v[2]);
  scs_float x = 0.0, y = 0.0, r;
  scs_int i;

  /* Check v membership in K_a */
  if (xh >= 0 && yh >= 0 &&
      POW_CONE_TOL + POWF(xh, a) * POWF(yh, (1 - a)) >= rh)
    return;

  /* Check -v membership in Polar K_a^* */
  if (xh <= 0 && yh <= 0 &&
      POW_CONE_TOL + POWF(-xh, a) * POWF(-yh, 1 - a) >=
          rh * POWF(a, a) * POWF(1 - a, 1 - a)) {
    v[0] = v[1] = v[2] = 0;
    return;
  }

  r = rh / 2;
  for (i = 0; i < POW_CONE_MAX_ITERS; ++i) {
    scs_float f, fp, dxdr, dydr;
    x = pow_calc_x(r, xh, rh, a);
    y = pow_calc_x(r, yh, rh, 1 - a);

    f = POWF(x, a) * POWF(y, (1 - a)) - r;
    if (ABS(f) < POW_CONE_TOL)
      break;

    dxdr = a * (rh - 2 * r) / (2 * x - xh);
    dydr = (1 - a) * (rh - 2 * r) / (2 * y - yh);
    fp = pow_calc_fp(x, y, dxdr, dydr, a);

    r = MAX(r - f / fp, 0);
    r = MIN(r, rh);
  }
  v[0] = x;
  v[1] = y;
  v[2] = (v[2] < 0) ? -r : r;
}

/* Project onto the primal K cone in the paper */
/* The r_y vector determines the INVERSE metric, ie, project under the
 * diag(r_y)^-1 norm.
 */
static scs_int proj_cone(scs_float *x, const ScsCone *k, ScsConeWork *c,
                         scs_int normalize, scs_float *r_y) {
  scs_int i, count = 0, status = 0;
  scs_float *r_box = SCS_NULL;
#ifdef USE_SPECTRAL_CONES
  SPECTRAL_TIMING(SCS(timer) spec_mat_proj_timer;)
#endif

  /* 1. Zero Cone */
  if (k->z) {
    memset(x, 0, k->z * sizeof(scs_float));
    count += k->z;
  }

  /* 2. Linear Cone (Non-negative orthant) */
  if (k->l) {
    for (i = count; i < count + k->l; ++i)
      x[i] = MAX(x[i], 0.0);
    count += k->l;
  }

  /* 3. Box Cone */
  if (k->bsize) {
    if (r_y)
      r_box = &(r_y[count]);
    c->box_t_warm_start = proj_box_cone(&(x[count]), k->bl, k->bu, k->bsize,
                                        c->box_t_warm_start, r_box);
    count += k->bsize;
  }

  /* 4. SOC */
  if (k->qsize) {
    for (i = 0; i < k->qsize; ++i) {
      proj_soc(&(x[count]), k->q[i]);
      count += k->q[i];
    }
  }

  /* 5. PSD Cone */
  if (k->ssize) {
    for (i = 0; i < k->ssize; ++i) {
#ifdef USE_SPECTRAL_CONES
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
#endif
      status = proj_semi_definite_cone(&(x[count]), k->s[i], c);
#ifdef USE_SPECTRAL_CONES
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
#endif
      if (status < 0)
        return status;
      count += get_sd_cone_size(k->s[i]);
    }
  }

  /* 6. Complex PSD Cone */
  if (k->cssize) {
    for (i = 0; i < k->cssize; ++i) {
      status = proj_complex_semi_definite_cone(&(x[count]), k->cs[i], c);
      if (status < 0)
        return status;
      count += get_csd_cone_size(k->cs[i]);
    }
  }

  if (k->ep || k->ed) { /* doesn't use r_y */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < k->ep + k->ed; ++i) {
      SCS(proj_pd_exp_cone)(&(x[count + 3 * i]), i < k->ep);
    }
    count += 3 * (k->ep + k->ed);
  }

  /* 8. Power Cone */
  if (k->psize) {
    scs_float v[3];
    scs_int idx;
    /* TODO: openmp not working well for power cone. */
    /*
    ifdef _OPENMP
    pragma omp parallel for private(v, idx)
    endif
    */
    for (i = 0; i < k->psize; ++i) { /* doesn't use r_y */
      idx = count + 3 * i;
      if (k->p[i] >= 0) {
        /* Primal power cone */
        proj_power_cone(&(x[idx]), k->p[i]);
      } else {
        /* Dual power cone via Moreau */
        v[0] = -x[idx];
        v[1] = -x[idx + 1];
        v[2] = -x[idx + 2];
        proj_power_cone(v, -k->p[i]);
        x[idx] += v[0];
        x[idx + 1] += v[1];
        x[idx + 2] += v[2];
      }
    }
    count += 3 * k->psize;
  }

#ifdef USE_SPECTRAL_CONES
  /* Spectral Cones (Logdet, Nuclear, Ell1, Sum Largest) */
  scs_int offset_log = 0;

  if (k->dsize) {
    for (i = 0; i < k->dsize; ++i) {
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
      status = SCS(proj_logdet_cone)(&(x[count]), k->d[i], c, offset_log,
                                     c->log_cone_warmstarts + i);
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
      if (status < 0)
        return status;
      offset_log += k->d[i] + 2;
      count += get_sd_cone_size(k->d[i]) + 2;
    }
  }
  if (k->nucsize) {
    for (i = 0; i < k->nucsize; ++i) {
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
      status = SCS(proj_nuclear_cone)(&(x[count]), k->nuc_m[i], k->nuc_n[i], c);
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
      if (status < 0)
        return status;
      count += k->nuc_m[i] * k->nuc_n[i] + 1;
    }
  }
  if (k->ell1_size) {
    for (i = 0; i < k->ell1_size; ++i) {
      SCS(proj_ell_one)(&(x[count]), k->ell1[i], c);
      count += k->ell1[i] + 1;
    }
  }
  if (k->sl_size) {
    for (i = 0; i < k->sl_size; ++i) {
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
      status =
          SCS(proj_sum_largest_evals)(&(x[count]), k->sl_n[i], k->sl_k[i], c);
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
      if (status < 0)
        return status;
      count += get_sd_cone_size(k->sl_n[i]) + 1;
    }
  }
#endif

  return 0;
}

/*
 * Public API
 */

ScsConeWork *SCS(init_cone)(ScsCone *k, scs_int m) {
  ScsConeWork *c = (ScsConeWork *)scs_calloc(1, sizeof(ScsConeWork));
  if (!c)
    return SCS_NULL;

  c->k = k;
  c->m = m;
  c->scaled_cones = 0;

  set_cone_boundaries(k, c);
  c->s = (scs_float *)scs_calloc(m, sizeof(scs_float));

  /* Set up workspaces if matrix cones are present */
  if ((k->ssize && k->s) || (k->cssize && k->cs)
#ifdef USE_SPECTRAL_CONES
      || (k->dsize) || (k->nucsize) || (k->sl_size)
#endif
  ) {
    if (set_up_cone_work_spaces(c, k) < 0) {
      SCS(finish_cone)(c);
      return SCS_NULL;
    }
  }

#ifdef USE_SPECTRAL_CONES
  if (k->ell1_size > 0 && k->ell1) {
    if (set_up_ell1_cone_work_space(c, k) < 0) {
      SCS(finish_cone)(c);
      return SCS_NULL;
    }
  }
#endif

  return c;
}

void scale_box_cone(ScsCone *k, ScsConeWork *c, const ScsScaling *scal) {
  if (k->bsize && k->bu && k->bl) {
    c->box_t_warm_start = 1.;
    if (scal) {
      /* also does some sanitizing */
      normalize_box_cone(k, &(scal->D[k->z + k->l]), k->bsize);
    }
  }
}

/* Outward facing cone projection routine, performs projection in-place.
   If normalize > 0 then will use normalized (equilibrated) cones if applicable.

   Moreau decomposition for R-norm projections:

    `x + R^{-1} \Pi_{C^*}^{R^{-1}} ( - R x ) = \Pi_C^R ( x )`

   where \Pi^R_C is the projection onto C under the R-norm:

    `||x||_R = \sqrt{x ' R x}`.

*/
scs_int SCS(proj_dual_cone)(scs_float *x, ScsConeWork *c,
                            const ScsScaling *scal, scs_float *r_y) {
  scs_int status, i;
  ScsCone *k = c->k;

  if (!c->scaled_cones) {
    /* Normalize box cone if applicable */
    if (k->bsize && k->bu && k->bl) {
      c->box_t_warm_start = 1.;
      if (scal)
        normalize_box_cone(k, &(scal->D[k->z + k->l]), k->bsize);
    }
    c->scaled_cones = 1;
  }

  /* Copy s = x */
  memcpy(c->s, x, c->m * sizeof(scs_float));

  /* x -> - Rx */
  for (i = 0; i < c->m; ++i) {
    x[i] *= r_y ? -r_y[i] : -1.0;
  }

  /* Project -x onto cone, x -> \Pi_{C^*}^{R^{-1}}(-x) under r_y metric */
  status = proj_cone(x, k, c, scal ? 1 : 0, r_y);

  /* Return x + R^{-1} \Pi_{C^*}^{R^{-1}} ( -x ) */
  for (i = 0; i < c->m; ++i) {
    if (r_y)
      x[i] = x[i] / r_y[i] + c->s[i];
    else
      x[i] += c->s[i];
  }

  return status;
}
