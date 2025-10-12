#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h" /* contains BLAS(X) macros and type info */
#include "util.h"
#if defined(_MSC_VER)
  /* MSVC: no C99 <complex.h> */
#else
  #include <complex.h>
#endif

// We need these definitions here to avoid including complex.h in scs_types.h.
// Including complex.h in scs_types.h causes issues when building with C++
// compilers. Reach out to Daniel Cederberg if you have any questions about the
// following 7 lines of code.
#if defined(_MSC_VER)
  /* MSVC C: use POD layout compatible with interleaved BLAS complex */
  typedef struct { double real, imag; } scs_blas_cdouble;
  typedef struct { float  real, imag; } scs_blas_cfloat;
  #ifndef SFLOAT
    #define SCS_BLAS_COMPLEX_TYPE scs_blas_cdouble
    #define SCS_BLAS_COMPLEX_CAST(x) ((scs_blas_cdouble *)(x))
  #else
    #define SCS_BLAS_COMPLEX_TYPE scs_blas_cfloat
    #define SCS_BLAS_COMPLEX_CAST(x) ((scs_blas_cfloat *)(x))
  #endif
#else
  /* GCC/Clang: keep using C99 _Complex */
  #ifndef SFLOAT
    #define SCS_BLAS_COMPLEX_CAST(x) ((double _Complex *)(x))
    #define SCS_BLAS_COMPLEX_TYPE double _Complex
  #else
    #define SCS_BLAS_COMPLEX_CAST(x) ((float _Complex *)(x))
    #define SCS_BLAS_COMPLEX_TYPE float _Complex
  #endif
#endif

#define BOX_CONE_MAX_ITERS (25)
#define POW_CONE_TOL (1e-9)
#define POW_CONE_MAX_ITERS (20)

/* Box cone limits (+ or -) taken to be INF */
#define MAX_BOX_VAL (1e15)

#ifdef USE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(syev)(const char *jobz, const char *uplo, blas_int *n, scs_float *a,
                blas_int *lda, scs_float *w, scs_float *work, blas_int *lwork,
                blas_int *info);
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

// Forward declare spectral matrix cone projections
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
#endif

/* Forward declare exponential cone projection routine.
 * Implemented in exp_cone.c.
 */
scs_float SCS(proj_pd_exp_cone)(scs_float *v0, scs_int primal);

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
  /* copy bu, bl */
  if (src->bsize > 1) {
    dest->bu = (scs_float *)scs_calloc(src->bsize - 1, sizeof(scs_float));
    memcpy(dest->bu, src->bu, (src->bsize - 1) * sizeof(scs_float));
    dest->bl = (scs_float *)scs_calloc(src->bsize - 1, sizeof(scs_float));
    memcpy(dest->bl, src->bl, (src->bsize - 1) * sizeof(scs_float));
  } else {
    dest->bu = SCS_NULL;
    dest->bl = SCS_NULL;
  }
  /* copy SOC */
  if (src->qsize > 0) {
    dest->q = (scs_int *)scs_calloc(src->qsize, sizeof(scs_int));
    memcpy(dest->q, src->q, src->qsize * sizeof(scs_int));
  } else {
    dest->q = SCS_NULL;
  }
  /* copy PSD cone */
  if (src->ssize > 0) {
    dest->s = (scs_int *)scs_calloc(src->ssize, sizeof(scs_int));
    memcpy(dest->s, src->s, src->ssize * sizeof(scs_int));
  } else {
    dest->s = SCS_NULL;
  }
  /* copy complex PSD cone */
  if (src->cssize > 0) {
    dest->cs = (scs_int *)scs_calloc(src->cssize, sizeof(scs_int));
    memcpy(dest->cs, src->cs, src->cssize * sizeof(scs_int));
  } else {
    dest->cs = SCS_NULL;
  }
  /* copy power cone */
  if (src->psize > 0) {
    dest->p = (scs_float *)scs_calloc(src->psize, sizeof(scs_float));
    memcpy(dest->p, src->p, src->psize * sizeof(scs_float));
  } else {
    dest->p = SCS_NULL;
  }
#ifdef USE_SPECTRAL_CONES
  /* copy logdet cone */
  if (src->dsize > 0) {
    dest->d = (scs_int *)scs_calloc(src->dsize, sizeof(scs_int));
    memcpy(dest->d, src->d, src->dsize * sizeof(scs_int));
  } else {
    dest->d = SCS_NULL;
  }
  /* copy nuclear norm cone*/
  if (src->nucsize > 0) {
    dest->nuc_m = (scs_int *)scs_calloc(src->nucsize, sizeof(scs_int));
    memcpy(dest->nuc_m, src->nuc_m, src->nucsize * sizeof(scs_int));
    dest->nuc_n = (scs_int *)scs_calloc(src->nucsize, sizeof(scs_int));
    memcpy(dest->nuc_n, src->nuc_n, src->nucsize * sizeof(scs_int));
  } else {
    dest->nuc_m = SCS_NULL;
    dest->nuc_n = SCS_NULL;
  }
  /* copy ell1-norm cone */
  if (src->ell1_size > 0) {
    dest->ell1 = (scs_int *)scs_calloc(src->ell1_size, sizeof(scs_int));
    memcpy(dest->ell1, src->ell1, src->ell1_size * sizeof(scs_int));
  } else {
    dest->ell1 = SCS_NULL;
  }
  /* copy sum-of-largest eigenvalues cone */
  if (src->sl_size > 0) {
    dest->sl_n = (scs_int *)scs_calloc(src->sl_size, sizeof(scs_int));
    memcpy(dest->sl_n, src->sl_n, src->sl_size * sizeof(scs_int));
    dest->sl_k = (scs_int *)scs_calloc(src->sl_size, sizeof(scs_int));
    memcpy(dest->sl_k, src->sl_k, src->sl_size * sizeof(scs_int));
  } else {
    dest->sl_n = SCS_NULL;
    dest->sl_k = SCS_NULL;
  }
#endif
}

/* set the vector of rho y terms, based on scale and cones */
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

  /* others */
  for (i = c->k->z; i < c->m; ++i) {
    r_y[i] = 1.0 / scale;
  }
}

/* the function f aggregates the entries within each cone */
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

static inline scs_int get_sd_cone_size(scs_int s) {
  return (s * (s + 1)) / 2;
}

static inline scs_int get_csd_cone_size(scs_int cs) {
  return cs * cs;
}

/*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size strictly
 * larger than 1, boundaries malloc-ed here so should be freed.
 */
void set_cone_boundaries(const ScsCone *k, ScsConeWork *c) {
  scs_int i, s_cone_sz, cs_cone_sz, count = 0;
#ifdef USE_SPECTRAL_CONES
  scs_int cone_boundaries_len = 1 + k->qsize + k->ssize + k->cssize + k->ed +
                                k->ep + k->psize + k->dsize + k->nucsize +
                                k->ell1_size + k->sl_size;
#else
  scs_int cone_boundaries_len =
      1 + k->qsize + k->ssize + k->cssize + k->ed + k->ep + k->psize;
#endif
  scs_int *b = (scs_int *)scs_calloc(cone_boundaries_len, sizeof(scs_int));
  /* cones that can be scaled independently */
  b[count] = k->z + k->l + k->bsize;
  count += 1; /* started at 0 now move to first entry */
  for (i = 0; i < k->qsize; ++i) {
    b[count + i] = k->q[i];
  }
  count += k->qsize;
  /* semidefinite cones */
  for (i = 0; i < k->ssize; ++i) {
    s_cone_sz = get_sd_cone_size(k->s[i]);
    b[count + i] = s_cone_sz;
  }
  count += k->ssize; /* size here not ssize * (ssize + 1) / 2 */
  for (i = 0; i < k->cssize; ++i) {
    cs_cone_sz = get_csd_cone_size(k->cs[i]);
    b[count + i] = cs_cone_sz;
  }
  count += k->cssize;
  /* exp cones */
  for (i = 0; i < k->ep + k->ed; ++i) {
    b[count + i] = 3;
  }
  count += k->ep + k->ed;
  /* power cones */
  for (i = 0; i < k->psize; ++i) {
    b[count + i] = 3;
  }
  count += k->psize;
#ifdef USE_SPECTRAL_CONES
  /* logdet cones */
  for (i = 0; i < k->dsize; ++i) {
    b[count + i] = get_sd_cone_size(k->d[i]) + 2;
  }
  count += k->dsize;
  /* nuclear norm cones */
  for (i = 0; i < k->nucsize; ++i) {
    b[count + i] = k->nuc_m[i] * k->nuc_n[i] + 1;
  }
  count += k->nucsize;
  /* ell1-norm cones */
  for (i = 0; i < k->ell1_size; ++i) {
    b[count + i] = k->ell1[i] + 1;
  }
  count += k->ell1_size;
  /* sum-of-largest eigenvalues cones */
  for (i = 0; i < k->sl_size; ++i) {
    b[count + i] = get_sd_cone_size(k->sl_n[i]) + 1;
  }
  count += k->sl_size;
#endif
  /* other cones */
  c->cone_boundaries = b;
  c->cone_boundaries_len = cone_boundaries_len;
}

static scs_int get_full_cone_dims(const ScsCone *k) {
  scs_int i, c = k->z + k->l + k->bsize;
  if (k->qsize) {
    for (i = 0; i < k->qsize; ++i) {
      c += k->q[i];
    }
  }
  if (k->ssize) {
    for (i = 0; i < k->ssize; ++i) {
      c += get_sd_cone_size(k->s[i]);
    }
  }
  if (k->cssize) {
    for (i = 0; i < k->cssize; ++i) {
      c += get_csd_cone_size(k->cs[i]);
    }
  }
  if (k->ed) {
    c += 3 * k->ed;
  }
  if (k->ep) {
    c += 3 * k->ep;
  }
  if (k->psize) {
    c += 3 * k->psize;
  }
#ifdef USE_SPECTRAL_CONES
  if (k->dsize) {
    for (i = 0; i < k->dsize; ++i) {
      c += get_sd_cone_size(k->d[i]) + 2;
    }
  }
  if (k->nucsize) {
    for (i = 0; i < k->nucsize; ++i) {
      c += k->nuc_m[i] * k->nuc_n[i] + 1;
    }
  }
  if (k->ell1_size) {
    for (i = 0; i < k->ell1_size; ++i) {
      c += k->ell1[i] + 1;
    }
  }
  if (k->sl_size) {
    for (i = 0; i < k->sl_size; ++i) {
      c += get_sd_cone_size(k->sl_n[i]) + 1;
    }
  }
#endif
  return c;
}

scs_int SCS(validate_cones)(const ScsData *d, const ScsCone *k) {
  scs_int i;
  if (get_full_cone_dims(k) != d->m) {
    scs_printf("cone dimensions %li not equal to num rows in A = m = %li\n",
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
#ifdef USE_LAPACK
  if (c->Xs) {
    scs_free(c->Xs);
  }
  if (c->cXs) {
    scs_free(c->cXs);
  }
  if (c->Z) {
    scs_free(c->Z);
  }
  if (c->cZ) {
    scs_free(c->cZ);
  }
  if (c->e) {
    scs_free(c->e);
  }
  if (c->isuppz) {
    scs_free(c->isuppz);
  }
  if (c->work) {
    scs_free(c->work);
  }
  if (c->iwork) {
    scs_free(c->iwork);
  }
  if (c->cwork) {
    scs_free(c->cwork);
  }
  if (c->rwork) {
    scs_free(c->rwork);
  }
#endif
  if (c->cone_boundaries) {
    scs_free(c->cone_boundaries);
  }
  if (c->s) {
    scs_free(c->s);
  }
#ifdef USE_SPECTRAL_CONES
  if (c->work_logdet) {
    scs_free(c->work_logdet);
  }
  if (c->saved_log_projs) {
    scs_free(c->saved_log_projs);
  }
  if (c->s_nuc) {
    scs_free(c->s_nuc);
  }
  if (c->u_nuc) {
    scs_free(c->u_nuc);
  }
  if (c->vt_nuc) {
    scs_free(c->vt_nuc);
  }
  if (c->work_nuc) {
    scs_free(c->work_nuc);
  }
  if (c->work_sum_of_largest) {
    scs_free(c->work_sum_of_largest);
  }
  if (c->log_cone_warmstarts) {
    scs_free(c->log_cone_warmstarts);
  }
  if (c->work_ell1) {
    scs_free(c->work_ell1);
  }
  if (c->work_ell1_proj) {
    scs_free(c->work_ell1_proj);
  }
#endif
  if (c) {
    scs_free(c);
  }
}

char *SCS(get_cone_header)(const ScsCone *k) {
  char *tmp = (char *)scs_malloc(512); /* sizeof(char) = 1 */
  scs_int i, soc_vars, sd_vars, csd_vars;
  sprintf(tmp, "cones: ");
  if (k->z) {
    sprintf(tmp + strlen(tmp), "\t  z: primal zero / dual free vars: %li\n",
            (long)k->z);
  }
  if (k->l) {
    sprintf(tmp + strlen(tmp), "\t  l: linear vars: %li\n", (long)k->l);
  }
  if (k->bsize) {
    sprintf(tmp + strlen(tmp), "\t  b: box cone vars: %li\n", (long)(k->bsize));
  }
  soc_vars = 0;
  if (k->qsize && k->q) {
    for (i = 0; i < k->qsize; i++) {
      soc_vars += k->q[i];
    }
    sprintf(tmp + strlen(tmp), "\t  q: soc vars: %li, qsize: %li\n",
            (long)soc_vars, (long)k->qsize);
  }
  sd_vars = 0;
  if (k->ssize && k->s) {
    for (i = 0; i < k->ssize; i++) {
      sd_vars += get_sd_cone_size(k->s[i]);
    }
    sprintf(tmp + strlen(tmp), "\t  s: psd vars: %li, ssize: %li\n",
            (long)sd_vars, (long)k->ssize);
  }
  csd_vars = 0;
  if (k->cssize && k->cs) {
    for (i = 0; i < k->cssize; i++) {
      csd_vars += get_csd_cone_size(k->cs[i]);
    }
    sprintf(tmp + strlen(tmp), "\t  cs: complex psd vars: %li, cssize: %li\n",
            (long)csd_vars, (long)k->cssize);
  }
  if (k->ep || k->ed) {
    sprintf(tmp + strlen(tmp), "\t  e: exp vars: %li, dual exp vars: %li\n",
            (long)(3 * k->ep), (long)(3 * k->ed));
  }
  if (k->psize && k->p) {
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

static scs_int set_up_sd_cone_work_space(ScsConeWork *c, const ScsCone *k) {
  scs_int i;
#ifdef USE_LAPACK
  blas_int n_max = 1;
  blas_int neg_one = -1;
  blas_int info = 0;
  scs_float wkopt = 0.0;
#if VERBOSITY > 0
#define _STR_EXPAND(tok) #tok
#define _STR(tok) _STR_EXPAND(tok)
  scs_printf("BLAS(func) = '%s'\n", _STR(BLAS(func)));
#endif

  // ----------------------------------------------------------------------
  //   compute max dimension needed for eigenvalue decomposition workspace
  //   (all cones appearing in the next section of code require eigenvalue
  //    decompositions)
  // ----------------------------------------------------------------------
  for (i = 0; i < k->ssize; ++i) {
    if (k->s[i] > n_max) {
      n_max = (blas_int)k->s[i];
    }
  }

#ifdef USE_SPECTRAL_CONES
  blas_int n_max_logdet = 1;
  blas_int n_logdet_total = 0;
  blas_int n_max_sl = 1;
  c->log_cone_warmstarts = (bool *)scs_calloc(k->dsize, sizeof(bool));

  for (i = 0; i < k->sl_size; ++i) {
    if (k->sl_n[i] > n_max_sl) {
      n_max_sl = (blas_int)k->sl_n[i];
    }
  }

  for (i = 0; i < k->dsize; ++i) {
    n_logdet_total += (blas_int)k->d[i];
    if (k->d[i] > n_max_logdet) {
      n_max_logdet = (blas_int)k->d[i];
    }
  }

  // --------------------------------------------------------------------------
  //         allocate workspace for logdeterminant cones
  // --------------------------------------------------------------------------
  if (k->dsize > 0) {
    c->work_logdet =
        (scs_float *)scs_calloc(22 * n_max_logdet + 122, sizeof(scs_float));
    c->saved_log_projs = (scs_float *)scs_calloc(2 * k->dsize + n_logdet_total,
                                                 sizeof(scs_float));

    if (!c->work_logdet || !c->saved_log_projs) {
      return -1;
    }
  }

  // ---------------------------------------------------------------
  //        allocate workspace for sum-of-largest-eigenvalues cone
  // ---------------------------------------------------------------
  if (k->sl_size > 0) {
    c->work_sum_of_largest =
        (scs_float *)scs_calloc(n_max_sl * n_max_sl, sizeof(scs_float));
    if (!c->work_sum_of_largest) {
      return -1;
    }
  }
  n_max = MAX(n_max, n_max_logdet);
  n_max = MAX(n_max, n_max_sl);
#endif
  // -----------------------------------------------------------------
  //         allocate eigenvector decomposition workspace
  // -----------------------------------------------------------------
  c->Xs = (scs_float *)scs_calloc(n_max * n_max, sizeof(scs_float));
  c->Z = (scs_float *)scs_calloc(n_max * n_max, sizeof(scs_float));
  c->e = (scs_float *)scs_calloc(n_max, sizeof(scs_float));

  /* workspace query */
  BLAS(syev)("Vectors", "Lower", &n_max, c->Xs, &n_max, SCS_NULL, &wkopt,
             &neg_one, &info);

  if (info != 0) {
    scs_printf("FATAL: syev workspace query failure, info = %li\n", (long)info);
    return -1;
  }

  c->lwork = (blas_int)(wkopt + 1); /* +1 for int casting safety */
  c->work = (scs_float *)scs_calloc(c->lwork, sizeof(scs_float));

  if (!c->Xs || !c->Z || !c->e || !c->work) {
    return -1;
  }

#ifdef USE_SPECTRAL_CONES
  // ------------------------------------------------------------------
  //              allocate memory for nuclear norm cone
  // ------------------------------------------------------------------
  if (k->nucsize > 0) {
    blas_int m_max_nuc = 1;
    blas_int n_max_nuc = 1;
    for (i = 0; i < k->nucsize; ++i) {
      if (k->nuc_m[i] > m_max_nuc) {
        m_max_nuc = k->nuc_m[i];
      }
      if (k->nuc_n[i] > n_max_nuc) {
        n_max_nuc = k->nuc_n[i];
      }
    }

    c->s_nuc = (scs_float *)scs_calloc(n_max_nuc, sizeof(scs_float));
    c->u_nuc =
        (scs_float *)scs_calloc(m_max_nuc * n_max_nuc, sizeof(scs_float));
    c->vt_nuc =
        (scs_float *)scs_calloc(n_max_nuc * n_max_nuc, sizeof(scs_float));

    if (!c->s_nuc || !c->u_nuc || !c->vt_nuc) {
      return -1;
    }

    // workspace query
    BLAS(gesvd)("S", "A", &m_max_nuc, &n_max_nuc, c->u_nuc, &m_max_nuc,
                c->s_nuc, c->u_nuc, &m_max_nuc, c->vt_nuc, &n_max_nuc, &wkopt,
                &neg_one, &info);

    c->lwork_nuc = (blas_int)(wkopt + 1); /* +1 for int casting safety */
    c->work_nuc = (scs_float *)scs_calloc(c->lwork_nuc, sizeof(scs_float));

    if (!c->lwork_nuc || !c->work_nuc) {
      return -1;
    }

    if (info != 0) {
      scs_printf("FATAL: gesvd workspace query failure, info = %li\n",
                 (long)info);
      return -1;
    }
  }
#endif

  return 0;
#else
  for (i = 0; i < k->ssize; i++) {
    if (k->s[i] > 1) {
      scs_printf(
          "FATAL: Cannot solve SDPs without linked blas+lapack libraries\n");
      scs_printf(
          "Install blas+lapack and re-compile SCS with blas+lapack library "
          "locations\n");
      return -1;
    }
  }
#ifdef USE_SPECTRAL_CONES
  bool has_spectral_cones = false;
  for (i = 0; i < k->dsize; i++) {
    if (k->d[i] > 1) {
      has_spectral_cones = true;
      break;
    }
  }

  for (i = 0; i < k->nucsize; i++) {
    if (k->nuc_m[i] > 1 || k->nuc_n[i] > 1) {
      has_spectral_cones = true;
      break;
    }
  }

  for (i = 0; i < k->sl_size; i++) {
    if (k->sl_n[i] > 1) {
      has_spectral_cones = true;
      break;
    }
  }

  if (has_spectral_cones) {
    scs_printf("FATAL: Cannot use spectral cones without linked blas+lapack "
               "libraries\n");
    scs_printf(
        "Install blas+lapack and re-compile SCS with blas+lapack library "
        "locations\n");
    return -1;
  }
#endif
  return 0;
#endif
}

/* size of X is get_sd_cone_size(n) */
static scs_int proj_semi_definite_cone(scs_float *X, const scs_int n,
                                       ScsConeWork *c) {
/* project onto the positive semi-definite cone */
#ifdef USE_LAPACK
  SCS(timer) _timer;
  scs_int i, first_idx;
  blas_int nb = (blas_int)n;
  blas_int ncols_z;
  blas_int nb_plus_one = (blas_int)(n + 1);
  blas_int one_int = 1;
  scs_float zero = 0., one = 1.;
  scs_float sqrt2 = SQRTF(2.0);
  scs_float sqrt2_inv = 1.0 / sqrt2;
  scs_float *Xs = c->Xs;
  scs_float *Z = c->Z;
  scs_float *e = c->e;
  scs_float *work = c->work;
  blas_int lwork = c->lwork;
  blas_int info = 0;
  scs_float sq_eig_pos;

#endif

  if (n == 0) {
    return 0;
  }
  if (n == 1) {
    X[0] = MAX(X[0], 0.);
    return 0;
  }

#ifdef USE_LAPACK

  /* copy lower triangular matrix into full matrix */
  for (i = 0; i < n; ++i) {
    memcpy(&(Xs[i * (n + 1)]), &(X[i * n - ((i - 1) * i) / 2]),
           (n - i) * sizeof(scs_float));
  }
  /*
     rescale so projection works, and matrix norm preserved
     see http://www.seas.ucla.edu/~vandenbe/publications/mlbook.pdf pg 3
   */
  /* scale diags by sqrt(2) */
  BLAS(scal)(&nb, &sqrt2, Xs, &nb_plus_one); /* not n_squared */

  /* Solve eigenproblem, reuse workspaces */
  BLAS(syev)("Vectors", "Lower", &nb, Xs, &nb, e, work, &lwork, &info);
  if (info != 0) {
    scs_printf("WARN: LAPACK syev error, info = %i\n", (int)info);
    if (info < 0) {
      return info;
    }
  }

  SCS(tic)(&_timer);

  first_idx = -1;
  /* e is eigvals in ascending order, find first entry > 0 */
  for (i = 0; i < n; ++i) {
    if (e[i] > 0) {
      first_idx = i;
      break;
    }
  }

  if (first_idx == -1) {
    /* there are no positive eigenvalues, set X to 0 and return */
    memset(X, 0, sizeof(scs_float) * get_sd_cone_size(n));
    return 0;
  }

  /* Z is matrix of eigenvectors with positive eigenvalues */
  memcpy(Z, &Xs[first_idx * n], sizeof(scs_float) * n * (n - first_idx));

  /* scale Z by sqrt(eig) */
  for (i = first_idx; i < n; ++i) {
    sq_eig_pos = SQRTF(e[i]);
    BLAS(scal)(&nb, &sq_eig_pos, &Z[(i - first_idx) * n], &one_int);
  }

  /* Xs = Z Z' = V E V' */
  ncols_z = (blas_int)(n - first_idx);
  BLAS(syrk)("Lower", "NoTrans", &nb, &ncols_z, &one, Z, &nb, &zero, Xs, &nb);

  /* undo rescaling: scale diags by 1/sqrt(2) */
  BLAS(scal)(&nb, &sqrt2_inv, Xs, &nb_plus_one); /* not n_squared */

  /* extract just lower triangular matrix */
  for (i = 0; i < n; ++i) {
    memcpy(&(X[i * n - ((i - 1) * i) / 2]), &(Xs[i * (n + 1)]),
           (n - i) * sizeof(scs_float));
  }
  return 0;

#else
  scs_printf("FAILURE: solving SDP but no blas/lapack libraries were found!\n");
  scs_printf("SCS will return nonsense!\n");
  SCS(scale_array)(X, NAN, get_sd_cone_size(n));
  return -1;
#endif
}

static scs_int set_up_csd_cone_work_space(ScsConeWork *c, const ScsCone *k) {
  scs_int i;
#ifdef USE_LAPACK
  blas_int n_max = 1;
  blas_int neg_one = -1;
  blas_int info = 0;
  scs_float abstol = -1.0;
  blas_int m = 0;
  scs_complex_float lcwork = {0.0};
  scs_float lrwork = 0.0;
  blas_int liwork = 0;
#if VERBOSITY > 0
#define _STR_EXPAND(tok) #tok
#define _STR(tok) _STR_EXPAND(tok)
  scs_printf("BLAS(func) = '%s'\n", _STR(BLASC(func)));
#endif
  /* eigenvector decomp workspace */
  for (i = 0; i < k->cssize; ++i) {
    if (k->cs[i] > n_max) {
      n_max = (blas_int)k->cs[i];
    }
  }
  c->cXs =
      (scs_complex_float *)scs_calloc(n_max * n_max, sizeof(scs_complex_float));
  c->cZ =
      (scs_complex_float *)scs_calloc(n_max * n_max, sizeof(scs_complex_float));
  c->e = (scs_float *)scs_calloc(n_max, sizeof(scs_float));
  c->isuppz = (blas_int *)scs_calloc(MAX(2, 2 * n_max), sizeof(blas_int));

  /* workspace query */
  BLASC(heevr)
  ("V", "A", "L", &n_max, SCS_BLAS_COMPLEX_CAST(c->cXs), &n_max, SCS_NULL,
   SCS_NULL, SCS_NULL, SCS_NULL, &abstol, &m, c->e,
   SCS_BLAS_COMPLEX_CAST(c->cZ), &n_max, c->isuppz,
   SCS_BLAS_COMPLEX_CAST(&lcwork), &neg_one, &lrwork, &neg_one, &liwork,
   &neg_one, &info);

  if (info != 0) {
    scs_printf("FATAL: heev workspace query failure, info = %li\n", (long)info);
    return -1;
  }
  c->lcwork = (blas_int)(lcwork[0]);
  c->lrwork = (blas_int)(lrwork);
  c->liwork = liwork;
  c->cwork =
      (scs_complex_float *)scs_calloc(c->lcwork, sizeof(scs_complex_float));
  c->rwork = (scs_float *)scs_calloc(c->lrwork, sizeof(scs_float));
  c->iwork = (blas_int *)scs_calloc(c->liwork, sizeof(blas_int));

  if (!c->cXs || !c->cZ || !c->e || !c->isuppz || !c->cwork || !c->rwork ||
      !c->iwork) {
    return -1;
  }
  return 0;
#else
  for (i = 0; i < k->cssize; i++) {
    if (k->cs[i] > 1) {
      scs_printf(
          "FATAL: Cannot solve SDPs without linked blas+lapack libraries\n");
      scs_printf(
          "Install blas+lapack and re-compile SCS with blas+lapack library "
          "locations\n");
      return -1;
    }
  }
  return 0;
#endif
}

/* size of X is get_csd_cone_size(n) */
static scs_int proj_complex_semi_definite_cone(scs_float *X, const scs_int n,
                                               ScsConeWork *c) {
/* project onto the positive semi-definite cone */
#ifdef USE_LAPACK
  scs_int i, first_idx;
  blas_int nb = (blas_int)n;
  blas_int ncols_z;
  blas_int nb_plus_one = (blas_int)(n + 1);
  blas_int one_int = 1;
  scs_float zero = 0., one = 1.;
  scs_complex_float csqrt2 = {0.0};
  csqrt2[0] = SQRTF(2.0);
  scs_complex_float csqrt2_inv = {0.0};
  csqrt2_inv[0] = 1.0 / csqrt2[0];
  ;
  scs_complex_float *cXs = c->cXs;
  scs_complex_float *cZ = c->cZ;
  scs_float *e = c->e;
  scs_float abstol = -1.0;
  blas_int m = 0;
  blas_int info = 0;
  scs_complex_float csq_eig_pos = {0.0};

#endif

  if (n == 0) {
    return 0;
  }
  if (n == 1) {
    X[0] = MAX(X[0], 0.);
    return 0;
  }

#ifdef USE_LAPACK

  /* copy lower triangular matrix into full matrix */
  for (i = 0; i < n - 1; ++i) {
    cXs[i * (n + 1)][0] = X[i * (2 * n - i)];
    cXs[i * (n + 1)][1] = 0.0;
    memcpy(&(cXs[i * (n + 1) + 1]), &(X[i * (2 * n - i) + 1]),
           2 * (n - i - 1) * sizeof(scs_float));
  }
  cXs[n * n - 1][0] = X[n * n - 1];
  cXs[n * n - 1][1] = 0.0;
  /*
     rescale so projection works, and matrix norm preserved
     see http://www.seas.ucla.edu/~vandenbe/publications/mlbook.pdf pg 3
   */
  /* scale diags by sqrt(2) */
  BLASC(scal)
  (&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2), SCS_BLAS_COMPLEX_CAST(cXs),
   &nb_plus_one); /* not n_squared */

  /* Solve eigenproblem, reuse workspaces */
  BLASC(heevr)
  ("V", "A", "L", &nb, SCS_BLAS_COMPLEX_CAST(cXs), &nb, SCS_NULL, SCS_NULL,
   SCS_NULL, SCS_NULL, &abstol, &m, e, SCS_BLAS_COMPLEX_CAST(cZ), &nb,
   c->isuppz, SCS_BLAS_COMPLEX_CAST(c->cwork), &c->lcwork, c->rwork, &c->lrwork,
   c->iwork, &c->liwork, &info);

  if (info != 0) {
    scs_printf("WARN: LAPACK heev error, info = %i\n", (int)info);
    if (info < 0) {
      return info;
    }
  }

  first_idx = -1;
  /* e is eigvals in ascending order, find first entry > 0 */
  for (i = 0; i < n; ++i) {
    if (e[i] > 0) {
      first_idx = i;
      break;
    }
  }

  if (first_idx == -1) {
    /* there are no positive eigenvalues, set X to 0 and return */
    memset(X, 0, sizeof(scs_float) * get_csd_cone_size(n));
    return 0;
  }

  /* cZ is matrix of all eigenvectors */
  /* scale cZ by sqrt(eig) */
  for (i = first_idx; i < n; ++i) {
    csq_eig_pos[0] = SQRTF(e[i]);
    BLASC(scal)
    (&nb, SCS_BLAS_COMPLEX_CAST(&csq_eig_pos),
     SCS_BLAS_COMPLEX_CAST(&cZ[i * n]), &one_int);
  }

  /* Xs = cZ cZ' = V E V' */
  ncols_z = (blas_int)(n - first_idx);
  BLASC(herk)
  ("Lower", "NoTrans", &nb, &ncols_z, &one,
   SCS_BLAS_COMPLEX_CAST(&cZ[first_idx * n]), &nb, &zero,
   SCS_BLAS_COMPLEX_CAST(cXs), &nb);

  /* undo rescaling: scale diags by 1/sqrt(2) */
  BLASC(scal)
  (&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2_inv), SCS_BLAS_COMPLEX_CAST(cXs),
   &nb_plus_one); /* not n_squared */

  /* extract just lower triangular matrix */
  for (i = 0; i < n - 1; ++i) {
    X[i * (2 * n - i)] = cXs[i * (n + 1)][0];
    memcpy(&(X[i * (2 * n - i) + 1]), &(cXs[i * (n + 1) + 1]),
           2 * (n - i - 1) * sizeof(scs_float));
  }
  X[n * n - 1] = cXs[n * n - 1][0];
  return 0;

#else
  scs_printf("FAILURE: solving SDP but no blas/lapack libraries were found!\n");
  scs_printf("SCS will return nonsense!\n");
  SCS(scale_array)(X, NAN, get_csd_cone_size(n));
  return -1;
#endif
}

static scs_float pow_calc_x(scs_float r, scs_float xh, scs_float rh,
                            scs_float a) {
  scs_float x = 0.5 * (xh + SQRTF(xh * xh + 4 * a * (rh - r) * r));
  return MAX(x, 1e-12);
}

static scs_float pow_calcdxdr(scs_float x, scs_float xh, scs_float rh,
                              scs_float r, scs_float a) {
  return a * (rh - 2 * r) / (2 * x - xh);
}

static scs_float pow_calc_f(scs_float x, scs_float y, scs_float r,
                            scs_float a) {
  return POWF(x, a) * POWF(y, (1 - a)) - r;
}

static scs_float pow_calc_fp(scs_float x, scs_float y, scs_float dxdr,
                             scs_float dydr, scs_float a) {
  return POWF(x, a) * POWF(y, (1 - a)) * (a * dxdr / x + (1 - a) * dydr / y) -
         1;
}

/*
 * Routine to scale the limits of the box cone by the scaling diagonal mat D > 0
 *
 *  want (t, s) \in K <==> (t', s') \in K'
 *
 *  (t', s') = (d0 * t, D s) (overloading D to mean D[1:])
 *    (up to scalar scaling factor which we can ignore due to conic property)
 *
 *   K = { (t, s) | t * l <= s <= t * u, t >= 0 } =>
 *       { (t, s) | d0 * t * D l / d0 <= D s <= d0 * t D u / d0, t >= 0 } =>
 *       { (t', s') | t' * l' <= s' <= t' u', t >= 0 } = K'
 *  where l' = D l  / d0, u' = D u / d0.
 */
static void normalize_box_cone(ScsCone *k, scs_float *D, scs_int bsize) {
  scs_int j;
  for (j = 0; j < bsize - 1; j++) {
    if (k->bu[j] >= MAX_BOX_VAL) {
      k->bu[j] = INFINITY;
    } else {
      k->bu[j] = D ? D[j + 1] * k->bu[j] / D[0] : k->bu[j];
    }
    if (k->bl[j] <= -MAX_BOX_VAL) {
      k->bl[j] = -INFINITY;
    } else {
      k->bl[j] = D ? D[j + 1] * k->bl[j] / D[0] : k->bl[j];
    }
  }
}

/* Project onto { (t, s) | t * l <= s <= t * u, t >= 0 }, Newton's method on t
   tx = [t; s], total length = bsize, under Euclidean metric 1/r_box.
*/
static scs_float proj_box_cone(scs_float *tx, const scs_float *bl,
                               const scs_float *bu, scs_int bsize,
                               scs_float t_warm_start, scs_float *r_box) {
  scs_float *x, gt, ht, t_prev, t = t_warm_start;
  scs_float rho_t = 1, *rho = SCS_NULL, r;
  scs_int iter, j;

  if (bsize == 1) { /* special case */
    tx[0] = MAX(tx[0], 0.0);
    return tx[0];
  }
  x = &(tx[1]);

  if (r_box) {
    rho_t = 1.0 / r_box[0];
    rho = &(r_box[1]);
  }

  /* should only require about 5 or so iterations, 1 or 2 if warm-started */
  for (iter = 0; iter < BOX_CONE_MAX_ITERS; iter++) {
    t_prev = t;
    gt = rho_t * (t - tx[0]); /* gradient */
    ht = rho_t;               /* hessian */
    for (j = 0; j < bsize - 1; j++) {
      r = rho ? 1.0 / rho[j] : 1.;
      if (x[j] > t * bu[j]) {
        gt += r * (t * bu[j] - x[j]) * bu[j]; /* gradient */
        ht += r * bu[j] * bu[j];              /* hessian */
      } else if (x[j] < t * bl[j]) {
        gt += r * (t * bl[j] - x[j]) * bl[j]; /* gradient */
        ht += r * bl[j] * bl[j];              /* hessian */
      }
    }
    t = MAX(t - gt / MAX(ht, 1e-8), 0.); /* newton step */
#if VERBOSITY > 3
    scs_printf("iter %i, t_new %1.3e, t_prev %1.3e, gt %1.3e, ht %1.3e\n", iter,
               t, t_prev, gt, ht);
    scs_printf("ABS(gt / (ht + 1e-6)) %.4e, ABS(t - t_prev) %.4e\n",
               ABS(gt / (ht + 1e-6)), ABS(t - t_prev));
#endif
    /* TODO: sometimes this check can fail (ie, declare convergence before it
     * should) if ht is very large, which can happen with some pathological
     * problems.
     */
    if (ABS(gt / MAX(ht, 1e-6)) < 1e-12 * MAX(t, 1.) ||
        ABS(t - t_prev) < 1e-11 * MAX(t, 1.)) {
      break;
    }
  }
  if (iter == BOX_CONE_MAX_ITERS) {
    scs_printf("warning: box cone proj hit maximum %i iters\n", (int)iter);
  }
  for (j = 0; j < bsize - 1; j++) {
    if (x[j] > t * bu[j]) {
      x[j] = t * bu[j];
    } else if (x[j] < t * bl[j]) {
      x[j] = t * bl[j];
    }
    /* x[j] unchanged otherwise */
  }
  tx[0] = t;

#if VERBOSITY > 3
  scs_printf("box cone iters %i\n", (int)iter + 1);
#endif
  return t;
}

/* project onto SOC of size q*/
static void proj_soc(scs_float *x, scs_int q) {
  if (q == 0) {
    return;
  }
  if (q == 1) {
    x[0] = MAX(x[0], 0.);
    return;
  }
  scs_float v1 = x[0];
  scs_float s = SCS(norm_2)(&(x[1]), q - 1);
  scs_float alpha = (s + v1) / 2.0;

  if (s <= v1) {
    return;
  } else if (s <= -v1) {
    memset(&(x[0]), 0, q * sizeof(scs_float));
  } else {
    x[0] = alpha;
    SCS(scale_array)(&(x[1]), alpha / s, q - 1);
  }
}

static void proj_power_cone(scs_float *v, scs_float a) {
  scs_float xh = v[0], yh = v[1], rh = ABS(v[2]);
  scs_float x = 0.0, y = 0.0, r;
  scs_int i;
  /* v in K_a */
  if (xh >= 0 && yh >= 0 &&
      POW_CONE_TOL + POWF(xh, a) * POWF(yh, (1 - a)) >= rh) {
    return;
  }

  /* -v in K_a^* */
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

    f = pow_calc_f(x, y, r, a);
    if (ABS(f) < POW_CONE_TOL) {
      break;
    }

    dxdr = pow_calcdxdr(x, xh, rh, r, a);
    dydr = pow_calcdxdr(y, yh, rh, r, (1 - a));
    fp = pow_calc_fp(x, y, dxdr, dydr, a);

    r = MAX(r - f / fp, 0);
    r = MIN(r, rh);
  }
  v[0] = x;
  v[1] = y;
  v[2] = (v[2] < 0) ? -(r) : (r);
}

/* project onto the primal K cone in the paper */
/* the r_y vector determines the INVERSE metric, ie, project under the
 * diag(r_y)^-1 norm.
 */
static scs_int proj_cone(scs_float *x, const ScsCone *k, ScsConeWork *c,
                         scs_int normalize, scs_float *r_y) {
  scs_int i, status;
  scs_int count = 0;
  scs_float *r_box = SCS_NULL;
#ifdef USE_SPECTRAL_CONES
  SPECTRAL_TIMING(SCS(timer) spec_mat_proj_timer;)
#endif

  if (k->z) { /* doesn't use r_y */
    /* project onto primal zero / dual free cone */
    memset(x, 0, k->z * sizeof(scs_float));
    count += k->z;
  }

  if (k->l) { /* doesn't use r_y */
    /* project onto positive orthant */
    for (i = count; i < count + k->l; ++i) {
      x[i] = MAX(x[i], 0.0);
    }
    count += k->l;
  }

  if (k->bsize) { /* DOES use r_y */
    if (r_y) {
      r_box = &(r_y[count]);
    }
    /* project onto box cone */
    c->box_t_warm_start = proj_box_cone(&(x[count]), k->bl, k->bu, k->bsize,
                                        c->box_t_warm_start, r_box);
    count += k->bsize; /* since b = (t,s), len(s) = bsize - 1 */
  }

  if (k->qsize && k->q) { /* doesn't use r_y */
    /* project onto second-order cones */
    for (i = 0; i < k->qsize; ++i) {
      proj_soc(&(x[count]), k->q[i]);
      count += k->q[i];
    }
  }

  if (k->ssize && k->s) { /* doesn't use r_y */
    /* project onto PSD cones */
    for (i = 0; i < k->ssize; ++i) {
#ifdef USE_SPECTRAL_CONES
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
#endif
      status = proj_semi_definite_cone(&(x[count]), k->s[i], c);
#ifdef USE_SPECTRAL_CONES
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
#endif
      if (status < 0) {
        return status;
      }
      count += get_sd_cone_size(k->s[i]);
    }
  }

  if (k->cssize && k->cs) { /* doesn't use r_y */
    /* project onto complex PSD cones */
    for (i = 0; i < k->cssize; ++i) {
      status = proj_complex_semi_definite_cone(&(x[count]), k->cs[i], c);
      if (status < 0) {
        return status;
      }
      count += get_csd_cone_size(k->cs[i]);
    }
  }

  if (k->ep || k->ed) { /* doesn't use r_y */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < k->ep + k->ed; ++i) {
      /* provided in exp_cone.c */
      SCS(proj_pd_exp_cone)(&(x[count + 3 * i]), i < k->ep);
    }
    count += 3 * (k->ep + k->ed);
  }
  if (k->psize && k->p) { /* doesn't use r_y */
    scs_float v[3];
    scs_int idx;
    /* don't use openmp for power cone
    ifdef _OPENMP
    pragma omp parallel for private(v, idx)
    endif
    */
    for (i = 0; i < k->psize; ++i) { /* doesn't use r_y */
      idx = count + 3 * i;
      if (k->p[i] >= 0) {
        /* primal power cone */
        proj_power_cone(&(x[idx]), k->p[i]);
      } else {
        /* dual power cone, using Moreau */
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
  scs_int offset_log_cone = 0; /* used for warmstarting log-cone projections */
  /* project onto logdet cones */
  if (k->dsize && k->d) {
    for (i = 0; i < k->dsize; ++i) {
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
      status = SCS(proj_logdet_cone)(&(x[count]), k->d[i], c, offset_log_cone,
                                     c->log_cone_warmstarts + i);
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
      offset_log_cone += k->d[i] + 2;
      if (status < 0) {
        return status;
      }
      count += (get_sd_cone_size(k->d[i]) + 2);
    }
  }
  /* project onto nuclear norm cones */
  if (k->nucsize && k->nuc_m && k->nuc_n) {
    for (i = 0; i < k->nucsize; ++i) {
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
      status = SCS(proj_nuclear_cone)(&(x[count]), k->nuc_m[i], k->nuc_n[i], c);
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
      if (status < 0) {
        return status;
      }
      count += (k->nuc_m[i] * k->nuc_n[i] + 1);
    }
  }
  /* project onto ell1-norm cones */
  if (k->ell1_size && k->ell1) {
    for (i = 0; i < k->ell1_size; ++i) {
      SCS(proj_ell_one)(&(x[count]), k->ell1[i], c);
      count += (k->ell1[i] + 1);
    }
  }
  /* project onto sum-of-largest eigenvalues cone */
  if (k->sl_size && k->sl_n && k->sl_k) {
    for (i = 0; i < k->sl_size; ++i) {
      SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
      status =
          SCS(proj_sum_largest_evals)(&(x[count]), k->sl_n[i], k->sl_k[i], c);
      SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                      SCS(tocq)(&spec_mat_proj_timer);)
      if (status < 0) {
        return status;
      }
      count += (get_sd_cone_size(k->sl_n[i]) + 1);
    }
  }

#endif
  /* project onto OTHER cones */

  return 0;
}

#ifdef USE_SPECTRAL_CONES
static scs_int set_up_ell1_cone_work_space(ScsConeWork *c, const ScsCone *k) {
  scs_int i;
  scs_int n_max = 0;

  if (k->ell1_size > 0) {
    for (i = 0; i < k->ell1_size; ++i) {
      n_max = MAX(k->ell1[i], n_max);
    }

    c->work_ell1 = (Value_index *)scs_calloc(n_max, sizeof(Value_index));
    c->work_ell1_proj = (scs_float *)scs_calloc(n_max + 1, sizeof(scs_float));
    if (!c->work_ell1 || !c->work_ell1_proj) {
      return -1;
    }
  }

  return 0;
}
#endif

ScsConeWork *SCS(init_cone)(ScsCone *k, scs_int m) {
  ScsConeWork *c = (ScsConeWork *)scs_calloc(1, sizeof(ScsConeWork));
  c->k = k;
  c->m = m;

  c->scaled_cones = 0;
  set_cone_boundaries(k, c);
  c->s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  if ((k->ssize && k->s)
#ifdef USE_SPECTRAL_CONES
      || (k->dsize && k->d) || (k->nucsize && k->nuc_m && k->nuc_n) ||
      (k->sl_size && k->sl_k && k->sl_n)
#endif
  ) {
    if (set_up_sd_cone_work_space(c, k) < 0) {
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

  if (k->cssize && k->cs) {
    if (set_up_csd_cone_work_space(c, k) < 0) {
      SCS(finish_cone)(c);
      return SCS_NULL;
    }
  }
  return c;
}

void scale_box_cone(ScsCone *k, ScsConeWork *c, ScsScaling *scal) {
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
scs_int SCS(proj_dual_cone)(scs_float *x, ScsConeWork *c, ScsScaling *scal,
                            scs_float *r_y) {
  scs_int status, i;
  ScsCone *k = c->k;

  if (!c->scaled_cones) {
    scale_box_cone(k, c, scal);
    c->scaled_cones = 1;
  }

  /* copy s = x */
  memcpy(c->s, x, c->m * sizeof(scs_float));

  /* x -> - Rx */
  for (i = 0; i < c->m; ++i) {
    x[i] *= r_y ? -r_y[i] : -1;
  }

  /* project -x onto cone, x -> \Pi_{C^*}^{R^{-1}}(-x) under r_y metric */
  status = proj_cone(x, k, c, scal ? 1 : 0, r_y);

  /* return x + R^{-1} \Pi_{C^*}^{R^{-1}} ( -x )  */
  for (i = 0; i < c->m; ++i) {
    if (r_y) {
      x[i] = x[i] / r_y[i] + c->s[i];
    } else {
      x[i] += c->s[i];
    }
  }

  return status;
}
