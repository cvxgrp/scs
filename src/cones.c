/*
 * Cone projection implementation.
 *
 * Handles projection onto the dual cone for all supported cone types:
 * zero/free, non-negative, box, second-order, semidefinite, complex
 * semidefinite, exponential, power, and (optionally) spectral cones.
 * Uses Moreau decomposition: proj_K*(x) = x + proj_K(-x).
 */

#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h" /* contains BLAS(X) macros and type info */
#include "util.h"

#include <float.h>
#include <stdio.h>
#include <string.h>
/* USE_PARALLEL_PSD is a build option, off by default, that compiles in the
 * parallel projection of the PSD blocks. It needs OpenMP and LAPACK. Every
 * other build, with OpenMP or not, projects the blocks one after another. It
 * is an option because it is only faster when the BLAS runs on one thread
 * inside the threads of SCS: see "PARALLEL PROJECTION OF THE PSD BLOCKS"
 * below for what it does and what it asks of the BLAS. */
#if defined(USE_PARALLEL_PSD) && defined(_OPENMP) && defined(USE_LAPACK)
#define PSD_PARALLEL (1)
#include <omp.h>
#include <stdlib.h> /* qsort */
#endif

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
/* Power cone projection: POW_CONE_TOL is the tolerance relative to |v| for
 * the Newton solve and the membership checks, POW_CONE_START_SLOPE bounds
 * |dx/dr| and |dy/dr| at r = rh for the Newton solve to start there, and
 * POW_CONE_EPS is the machine epsilon. All must suit the working precision;
 * in particular the first Newton step from rh is only representable when
 * POW_CONE_START_SLOPE^2 * POW_CONE_EPS is small, hence the float value. */
#ifdef SFLOAT
#define POW_CONE_TOL (1e-6)
#define POW_CONE_START_SLOPE (10)
#define POW_CONE_EPS (FLT_EPSILON)
#else
#define POW_CONE_TOL (1e-12)
#define POW_CONE_START_SLOPE (1e3)
#define POW_CONE_EPS (DBL_EPSILON)
#endif
#define POW_CONE_MAX_ITERS (50)

/* Box cone limits (+ or -) taken to be INF */
#define MAX_BOX_VAL (1e15)

#ifdef USE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif

/* LAPACK / BLAS Function Prototypes */
/* The eigendecompositions of the PSD projection: divide and conquer, with the
 * eigenvectors returned in place of the matrix (see proj_semi_definite_cone) */
void BLAS(syevd)(const char *jobz, const char *uplo, blas_int *n, scs_float *a,
                 blas_int *lda, scs_float *w, scs_float *work, blas_int *lwork,
                 blas_int *iwork, blas_int *liwork, blas_int *info);

void BLASC(heevd)(const char *jobz, const char *uplo, blas_int *n,
                  SCS_BLAS_COMPLEX_TYPE *a, blas_int *lda, scs_float *w,
                  SCS_BLAS_COMPLEX_TYPE *cwork, blas_int *lcwork,
                  scs_float *rwork, blas_int *lrwork, blas_int *iwork,
                  blas_int *liwork, blas_int *info);

void BLAS(syrk)(const char *uplo, const char *trans, const blas_int *n,
                const blas_int *k, const scs_float *alpha, const scs_float *a,
                const blas_int *lda, const scs_float *beta, scs_float *c,
                const blas_int *ldc);

void BLASC(herk)(const char *uplo, const char *trans, const blas_int *n,
                 const blas_int *k, const scs_float *alpha,
                 const SCS_BLAS_COMPLEX_TYPE *a, const blas_int *lda,
                 const scs_float *beta, SCS_BLAS_COMPLEX_TYPE *c,
                 const blas_int *ldc);

void BLAS(scal)(const blas_int *n, const scs_float *sa, scs_float *sx,
                const blas_int *incx);
void BLASC(scal)(const blas_int *n, const SCS_BLAS_COMPLEX_TYPE *sa,
                 SCS_BLAS_COMPLEX_TYPE *sx, const blas_int *incx);

/* Fallbacks for when syevd / heevd report no convergence, see
 * proj_semi_definite_cone. syev is also used by the spectral cones. */
void BLAS(syev)(const char *jobz, const char *uplo, blas_int *n, scs_float *a,
                blas_int *lda, scs_float *w, scs_float *work, blas_int *lwork,
                blas_int *info);

void BLASC(heev)(const char *jobz, const char *uplo, blas_int *n,
                 SCS_BLAS_COMPLEX_TYPE *a, blas_int *lda, scs_float *w,
                 SCS_BLAS_COMPLEX_TYPE *work, blas_int *lwork, scs_float *rwork,
                 blas_int *info);

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

/* ======================== Memory Management ======================== */

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

scs_int SCS(deep_copy_cone)(ScsCone *dest, const ScsCone *src) {
  memset(dest, 0, sizeof(*dest));
  dest->z = src->z;
  dest->l = src->l;
  dest->bsize = src->bsize;
  dest->qsize = src->qsize;
  dest->ssize = src->ssize;
  dest->cssize = src->cssize;
  dest->ep = src->ep;
  dest->ed = src->ed;
  dest->psize = src->psize;
#ifdef USE_SPECTRAL_CONES
  dest->dsize = src->dsize;
  dest->nucsize = src->nucsize;
  dest->ell1_size = src->ell1_size;
  dest->sl_size = src->sl_size;
#endif

  /* Box cone */
  if (src->bsize > 1) {
    dest->bu = (scs_float *)scs_calloc(src->bsize - 1, sizeof(scs_float));
    dest->bl = (scs_float *)scs_calloc(src->bsize - 1, sizeof(scs_float));
    if (!dest->bu || !dest->bl) {
      return 0;
    }
    memcpy(dest->bu, src->bu, (src->bsize - 1) * sizeof(scs_float));
    memcpy(dest->bl, src->bl, (src->bsize - 1) * sizeof(scs_float));
  } else {
    dest->bu = SCS_NULL;
    dest->bl = SCS_NULL;
  }

  /* SOC */
  if (src->qsize > 0) {
    dest->q = (scs_int *)scs_calloc(src->qsize, sizeof(scs_int));
    if (!dest->q) {
      return 0;
    }
    memcpy(dest->q, src->q, src->qsize * sizeof(scs_int));
  } else {
    dest->q = SCS_NULL;
  }

  /* PSD */
  if (src->ssize > 0) {
    dest->s = (scs_int *)scs_calloc(src->ssize, sizeof(scs_int));
    if (!dest->s) {
      return 0;
    }
    memcpy(dest->s, src->s, src->ssize * sizeof(scs_int));
  } else {
    dest->s = SCS_NULL;
  }

  /* Complex PSD */
  if (src->cssize > 0) {
    dest->cs = (scs_int *)scs_calloc(src->cssize, sizeof(scs_int));
    if (!dest->cs) {
      return 0;
    }
    memcpy(dest->cs, src->cs, src->cssize * sizeof(scs_int));
  } else {
    dest->cs = SCS_NULL;
  }

  /* Power */
  if (src->psize > 0) {
    dest->p = (scs_float *)scs_calloc(src->psize, sizeof(scs_float));
    if (!dest->p) {
      return 0;
    }
    memcpy(dest->p, src->p, src->psize * sizeof(scs_float));
  } else {
    dest->p = SCS_NULL;
  }

#ifdef USE_SPECTRAL_CONES
  /* Logdet */
  if (src->dsize > 0) {
    dest->d = (scs_int *)scs_calloc(src->dsize, sizeof(scs_int));
    if (!dest->d) {
      return 0;
    }
    memcpy(dest->d, src->d, src->dsize * sizeof(scs_int));
  } else {
    dest->d = SCS_NULL;
  }

  /* Nuclear */
  if (src->nucsize > 0) {
    dest->nuc_m = (scs_int *)scs_calloc(src->nucsize, sizeof(scs_int));
    dest->nuc_n = (scs_int *)scs_calloc(src->nucsize, sizeof(scs_int));
    if (!dest->nuc_m || !dest->nuc_n) {
      return 0;
    }
    memcpy(dest->nuc_m, src->nuc_m, src->nucsize * sizeof(scs_int));
    memcpy(dest->nuc_n, src->nuc_n, src->nucsize * sizeof(scs_int));
  } else {
    dest->nuc_m = SCS_NULL;
    dest->nuc_n = SCS_NULL;
  }

  /* Ell1 */
  if (src->ell1_size > 0) {
    dest->ell1 = (scs_int *)scs_calloc(src->ell1_size, sizeof(scs_int));
    if (!dest->ell1) {
      return 0;
    }
    memcpy(dest->ell1, src->ell1, src->ell1_size * sizeof(scs_int));
  } else {
    dest->ell1 = SCS_NULL;
  }

  /* Sum Largest Eigenvalues */
  if (src->sl_size > 0) {
    dest->sl_n = (scs_int *)scs_calloc(src->sl_size, sizeof(scs_int));
    dest->sl_k = (scs_int *)scs_calloc(src->sl_size, sizeof(scs_int));
    if (!dest->sl_n || !dest->sl_k) {
      return 0;
    }
    memcpy(dest->sl_n, src->sl_n, src->sl_size * sizeof(scs_int));
    memcpy(dest->sl_k, src->sl_k, src->sl_size * sizeof(scs_int));
  } else {
    dest->sl_n = SCS_NULL;
    dest->sl_k = SCS_NULL;
  }
#endif
  return 1;
}

void SCS(finish_cone)(ScsConeWork *c) {
  if (!c)
    return;
#ifdef USE_LAPACK
  if (c->psd_tw) {
    scs_int t, j;
    for (t = 0; t < MAX(c->psd_nthreads, 1); ++t) {
      ScsPsdWork *w = &(c->psd_tw[t]);
      void *bufs[] = {w->Xs,  w->Z,     w->e,    w->work,
                      w->gamma, w->cXs, w->cwork, w->iwork};
      for (j = 0; j < (scs_int)(sizeof(bufs) / sizeof(bufs[0])); ++j) {
        if (bufs[j])
          scs_free(bufs[j]);
      }
    }
    scs_free(c->psd_tw);
  }
  if (c->psd_lsz)
    scs_free(c->psd_lsz);
#endif
  if (c->psd_blk)
    scs_free(c->psd_blk);
  if (c->psd_tptr)
    scs_free(c->psd_tptr);
  if (c->psd_off)
    scs_free(c->psd_off);
  if (c->cone_boundaries)
    scs_free(c->cone_boundaries);
  if (c->cone_boundaries_psd_n)
    scs_free(c->cone_boundaries_psd_n);
  if (c->r_box_inv)
    scs_free(c->r_box_inv);
  if (c->box_kinks)
    scs_free(c->box_kinks);
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

/* ========================= Cone Utilities ============================ */

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
    r_y[i] = 1.0 / (ZERO_CONE_RY_FACTOR * scale);
  }
  /* Remaining cones */
  for (i = c->k->z; i < c->m; ++i) {
    r_y[i] = 1.0 / scale;
  }
  /* Box cone t entry: weighted like an equality row (BOX_T_RY_FACTOR in
   * glbopts.h says why). A smaller rho_y on this entry acts in the two places
   * rho_y enters:
   * - The linear system. In the usual form the t row of A is zero, so the
   *   row is decoupled from x: it gives y_t = v_t - tau b_t / rho_y, and adds
   *   b_t^2 / rho_y to the leading coefficient of the quadratic for tau
   *   (root_plus), as a zero-cone row with b != 0 does.
   * - The box projection, whose metric weighs t by rho_t = 1 / rho_y
   *   (proj_cone hands proj_box_cone the slice of r_y that starts at this
   *   entry). project_cones hands proj_dual_cone 2 y_t - v_t, which it
   *   negates and multiplies by rho_y, so with y_t from the linear system as
   *   above (tau its root_plus value) the t input is
   *   t0 = -rho_y (2 y_t - v_t) = 2 tau b_t - rho_y v_t. At t* = K / H (see
   *   proj_box_cone) the clipped rows pull t away from t0 by
   *   sum_j r_j b_j (x_j - t0 b_j) / H, with H = rho_t + S and
   *   S = sum_j r_j b_j^2 over them (b_j the bound row j is clipped to). A
   *   rho_t BOX_T_RY_FACTOR times larger shrinks that pull by
   *   (rho_t + S) / (BOX_T_RY_FACTOR rho_t + S), by up to BOX_T_RY_FACTOR when
   *   S does not dominate, so the box rows are clipped to a t nearer tau b_t.
   * Everything else reads the entry from diag_r: the linear system,
   * root_plus, the remap of v at a scale update and warm starts. With
   * adaptive_diag_scale the row's own multiplier divides it, as it divides a
   * zero-cone row's. So this is the one place it is set.
   * A box cone of size 1 is t >= 0 alone, with no bounds for t to scale: a
   * nonnegative entry, which keeps 1 / scale. */
  if (c->k->bsize > 1) {
    r_y[c->k->z + c->k->l] = 1.0 / (BOX_T_RY_FACTOR * scale);
  }
}

/* svec packed index of entry (i, j), i >= j, in an n x n lower
 * triangle stored column-major. */
static scs_int svec_idx(scs_int i, scs_int j, scs_int n) {
  return j * n - (j * (j - 1)) / 2 + (i - j);
}

/* Fit the rank-one form w_ij = delta_i * delta_j to a positive profile
 * on one PSD block (svec packed, length n(n+1)/2) by least squares in
 * log space, and overwrite the block with the fitted values. This is
 * the richest svec-diagonal metric that maps the PSD cone to itself
 * (the diagonal congruence X -> diag(delta) X diag(delta)), so unlike
 * generic cones the block need not collapse to a single scalar.
 *
 * Normal equations: minimizing sum_{i<=j} (l_i + l_j - y_ij)^2 gives
 * A'A = (n+2) I + 11', inverted in closed form by Sherman-Morrison.
 * The fitted delta_i are clamped so every product w_ij stays within
 * [DIAG_SCALE_MULT_MIN, DIAG_SCALE_MULT_MAX], preserving rank-one
 * structure (clamping products directly would break it). */
static void fit_psd_block_weights(scs_float *vec, scs_int n) {
  scs_int i, j;
  scs_float alpha = (scs_float)(n + 2);
  scs_float rhs_sum = 0., corr, li;
  scs_float lo = SQRTF(DIAG_SCALE_MULT_MIN), hi = SQRTF(DIAG_SCALE_MULT_MAX);
  scs_float *ell = (scs_float *)scs_calloc(n, sizeof(scs_float));
  if (!ell) {
    return; /* leave the block as it is */
  }
  /* Floor entries relative to the block max: (near-)zero rows carry no
   * information about the block's scale structure (e.g. unconstrained
   * off-diagonal svec rows are exactly zero in the equilibration
   * profile) and an absolute floor lets them dominate the fit. */
  {
    scs_float vmax = 0., vfloor;
    scs_int sz = (n * (n + 1)) / 2;
    for (i = 0; i < sz; ++i) {
      vmax = MAX(vmax, vec[i]);
    }
    vfloor = MAX(1e-3 * vmax, 1e-12);
    for (i = 0; i < sz; ++i) {
      vec[i] = MAX(vec[i], vfloor);
    }
  }

  /* rhs_k = 2 y_kk + sum_{j != k} y_kj */
  for (j = 0; j < n; ++j) {
    for (i = j; i < n; ++i) {
      scs_float y = log(vec[svec_idx(i, j, n)]);
      if (i == j) {
        ell[i] += 2. * y;
      } else {
        ell[i] += y;
        ell[j] += y;
      }
    }
  }
  for (i = 0; i < n; ++i) {
    rhs_sum += ell[i];
  }
  /* (alpha I + 11')^{-1} rhs = rhs/alpha - 1 (1'rhs) / (alpha(alpha+n)) */
  corr = rhs_sum / (alpha * (alpha + (scs_float)n));
  for (i = 0; i < n; ++i) {
    li = exp(ell[i] / alpha - corr);
    ell[i] = MIN(MAX(li, lo), hi);
  }
  for (j = 0; j < n; ++j) {
    for (i = j; i < n; ++i) {
      vec[svec_idx(i, j, n)] = ell[i] * ell[j];
    }
  }
  scs_free(ell);
}

/* The function f aggregates the entries within each cone. With
 * psd_rank1 set, real PSD blocks instead get the rank-one (diagonal
 * congruence) fit, which is cone-invariant and strictly more
 * expressive than a single scalar. The fit is only safe for the
 * one-shot dynamic metric update: applied inside the iterated Ruiz
 * equilibration loop it compounds pass over pass instead of
 * contracting and destroys the scaling (mcp100 diverges), so the
 * equilibration call sites keep the scalar collapse. */
void SCS(enforce_cone_boundaries)(const ScsConeWork *c, scs_float *vec,
                                  scs_float (*f)(const scs_float *, scs_int),
                                  scs_int psd_rank1) {
  scs_int i, j, delta;
  scs_int count = c->cone_boundaries[0];
  scs_float wrk;
  for (i = 1; i < c->cone_boundaries_len; ++i) {
    delta = c->cone_boundaries[i];
    if (psd_rank1 && c->cone_boundaries_psd_n &&
        c->cone_boundaries_psd_n[i] > 1) {
      fit_psd_block_weights(&(vec[count]), c->cone_boundaries_psd_n[i]);
    } else {
      wrk = f(&(vec[count]), delta);
      for (j = count; j < count + delta; ++j) {
        vec[j] = wrk;
      }
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
  scs_int *psd_n = (scs_int *)scs_calloc(total_cones + 1, sizeof(scs_int));

  /* Cones that can be scaled independently */
  b[count++] = k->z + k->l + k->bsize;
  for (i = 0; i < k->qsize; ++i)
    b[count++] = k->q[i];
  for (i = 0; i < k->ssize; ++i) {
    psd_n[count] = k->s[i];
    b[count++] = get_sd_cone_size(k->s[i]);
  }
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
  c->cone_boundaries_psd_n = psd_n;
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

/* If buf is NULL, only accumulate the length that would be written.
 * Otherwise append chunk, including its trailing '\0', at the current offset.
 */
static void append_to_header(char *buf, size_t *len, const char *chunk) {
  size_t chunk_len = strlen(chunk);
  if (buf) {
    memcpy(buf + *len, chunk, chunk_len + 1);
  }
  *len += chunk_len;
}

/* Format the cone summary in one shared code path.
 *
 * When buf is NULL, this function performs a sizing pass only and updates
 * *len with the number of characters required for the final string.
 * When buf is non-NULL, it writes the same formatted output into buf and
 * advances *len as it appends each line.
 *
 * This supports a two-pass implementation in get_cone_header:
 * 1. call with buf == NULL to compute the exact allocation size
 * 2. allocate once
 * 3. call again with buf != NULL to render the final string
 */
static void format_cone_header(const ScsCone *k, char *buf, size_t *len) {
  char line[128];
  scs_int i, count;
#ifdef USE_SPECTRAL_CONES
  scs_int ell1_vars, log_vars, nuc_vars, sl_vars;
#endif

  append_to_header(buf, len, "cones: ");
  if (k->z) {
    sprintf(line, "\t  z: primal zero / dual free vars: %li\n", (long)k->z);
    append_to_header(buf, len, line);
  }
  if (k->l) {
    sprintf(line, "\t  l: linear vars: %li\n", (long)k->l);
    append_to_header(buf, len, line);
  }
  if (k->bsize) {
    sprintf(line, "\t  b: box cone vars: %li\n", (long)k->bsize);
    append_to_header(buf, len, line);
  }
  if (k->qsize) {
    count = 0;
    for (i = 0; i < k->qsize; ++i)
      count += k->q[i];
    sprintf(line, "\t  q: soc vars: %li, qsize: %li\n", (long)count,
            (long)k->qsize);
    append_to_header(buf, len, line);
  }
  if (k->ssize) {
    count = 0;
    for (i = 0; i < k->ssize; ++i)
      count += get_sd_cone_size(k->s[i]);
    sprintf(line, "\t  s: psd vars: %li, ssize: %li\n", (long)count,
            (long)k->ssize);
    append_to_header(buf, len, line);
  }
  if (k->cssize) {
    count = 0;
    for (i = 0; i < k->cssize; ++i)
      count += get_csd_cone_size(k->cs[i]);
    sprintf(line, "\t  cs: complex psd vars: %li, cssize: %li\n", (long)count,
            (long)k->cssize);
    append_to_header(buf, len, line);
  }
  if (k->ep || k->ed) {
    sprintf(line, "\t  e: exp vars: %li, dual exp vars: %li\n",
            (long)(3 * k->ep), (long)(3 * k->ed));
    append_to_header(buf, len, line);
  }
  if (k->psize) {
    sprintf(line, "\t  p: primal + dual power vars: %li\n",
            (long)(3 * k->psize));
    append_to_header(buf, len, line);
  }
#ifdef USE_SPECTRAL_CONES
  log_vars = 0;
  if (k->dsize && k->d) {
    for (i = 0; i < k->dsize; i++) {
      log_vars += get_sd_cone_size(k->d[i]) + 2;
    }
    sprintf(line, "\t  d: logdet vars: %li, dsize: %li\n", (long)log_vars,
            (long)k->dsize);
    append_to_header(buf, len, line);
  }
  nuc_vars = 0;
  if (k->nucsize && k->nuc_m && k->nuc_n) {
    for (i = 0; i < k->nucsize; i++) {
      nuc_vars += k->nuc_m[i] * k->nuc_n[i] + 1;
    }
    sprintf(line, "\t  nuc: nuclear vars: %li, nucsize: %li\n",
            (long)nuc_vars, (long)k->nucsize);
    append_to_header(buf, len, line);
  }
  ell1_vars = 0;
  if (k->ell1_size && k->ell1) {
    for (i = 0; i < k->ell1_size; ++i) {
      ell1_vars += k->ell1[i] + 1;
    }
    sprintf(line, "\t  ell1: ell1 vars: %li, ell1_size: %li\n",
            (long)ell1_vars, (long)k->ell1_size);
    append_to_header(buf, len, line);
  }
  sl_vars = 0;
  if (k->sl_size && k->sl_n) {
    for (i = 0; i < k->sl_size; ++i) {
      sl_vars += get_sd_cone_size(k->sl_n[i]) + 1;
    }
    sprintf(line, "\t  sl: sl vars: %li, sl_size: %li\n", (long)sl_vars,
            (long)k->sl_size);
    append_to_header(buf, len, line);
  }
#endif
}

char *SCS(get_cone_header)(const ScsCone *k) {
  char *tmp;
  size_t len = 0;

  format_cone_header(k, SCS_NULL, &len);

  tmp = (char *)scs_malloc(len + 1);
  if (!tmp) {
    return SCS_NULL;
  }

  len = 0;
  format_cone_header(k, tmp, &len);
  return tmp;
}

/* ========================== Validation =============================== */

scs_int SCS(validate_cones)(const ScsData *d, const ScsCone *k) {
  scs_int i;
  scs_int cone_dims;
  if (k->z < 0) {
    scs_printf("free cone dimension error\n");
    return -1;
  }
  if (k->l < 0) {
    scs_printf("lp cone dimension error\n");
    return -1;
  }
  if (k->bsize < 0) {
    scs_printf("box cone dimension error\n");
    return -1;
  }
  if (k->bsize > 1) {
    if (!k->bl || !k->bu) {
      scs_printf("box cone bounds missing\n");
      return -1;
    }
    for (i = 0; i < k->bsize - 1; ++i) {
      if (isnan(k->bl[i]) || isnan(k->bu[i])) {
        scs_printf("box cone bounds cannot be NaN\n");
        return -1;
      }
      if (k->bl[i] == INFINITY || k->bu[i] == -INFINITY) {
        scs_printf("box cone bounds use invalid infinity direction\n");
        return -1;
      }
      if (k->bl[i] > k->bu[i]) {
        scs_printf("infeasible: box lower bound larger than upper bound\n");
        return -1;
      }
    }
  }
  if (k->qsize < 0) {
    scs_printf("soc cone dimension error\n");
    return -1;
  }
  if (k->qsize > 0) {
    if (!k->q) {
      scs_printf("soc cone array missing\n");
      return -1;
    }
    for (i = 0; i < k->qsize; ++i) {
      if (k->q[i] < 0) {
        scs_printf("soc cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->ssize < 0) {
    scs_printf("sd cone dimension error\n");
    return -1;
  }
  if (k->ssize > 0) {
    if (!k->s) {
      scs_printf("sd cone array missing\n");
      return -1;
    }
    for (i = 0; i < k->ssize; ++i) {
      if (k->s[i] < 0) {
        scs_printf("sd cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->cssize < 0) {
    scs_printf("complex psd cone dimension error\n");
    return -1;
  }
  if (k->cssize > 0) {
    if (!k->cs) {
      scs_printf("complex psd cone array missing\n");
      return -1;
    }
    for (i = 0; i < k->cssize; ++i) {
      if (k->cs[i] < 0) {
        scs_printf("complex psd cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->ed < 0) {
    scs_printf("ed cone dimension error\n");
    return -1;
  }
  if (k->ep < 0) {
    scs_printf("ep cone dimension error\n");
    return -1;
  }
  if (k->psize < 0) {
    scs_printf("power cone dimension error\n");
    return -1;
  }
  if (k->psize > 0) {
    if (!k->p) {
      scs_printf("power cone array missing\n");
      return -1;
    }
    for (i = 0; i < k->psize; ++i) {
      if (!isfinite(k->p[i]) || k->p[i] < -1 || k->p[i] > 1) {
        scs_printf("power cone error, values must be finite and in [-1,1]\n");
        return -1;
      }
    }
  }
#ifdef USE_SPECTRAL_CONES
  if (k->dsize < 0) {
    scs_printf("logdet cone dimension error\n");
    return -1;
  }
  if (k->dsize > 0) {
    if (!k->d) {
      scs_printf("logdet cone array missing\n");
      return -1;
    }
    for (i = 0; i < k->dsize; ++i) {
      if (k->d[i] < 1) {
        scs_printf("logdet cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->nucsize < 0) {
    scs_printf("nuclear cone dimension error\n");
    return -1;
  }
  if (k->nucsize > 0) {
    if (!k->nuc_m || !k->nuc_n) {
      scs_printf("nuclear cone arrays missing\n");
      return -1;
    }
    for (i = 0; i < k->nucsize; ++i) {
      if (k->nuc_m[i] < 1 || k->nuc_n[i] < 1 || k->nuc_n[i] > k->nuc_m[i]) {
        scs_printf("nuclear norm cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->ell1_size < 0) {
    scs_printf("ell1 cone dimension error\n");
    return -1;
  }
  if (k->ell1_size > 0) {
    if (!k->ell1) {
      scs_printf("ell1 cone array missing\n");
      return -1;
    }
    for (i = 0; i < k->ell1_size; ++i) {
      if (k->ell1[i] < 1) {
        scs_printf("ell1 cone dimension error\n");
        return -1;
      }
    }
  }
  if (k->sl_size < 0) {
    scs_printf("sum-of-largest-eigenvalues cone dimension error\n");
    return -1;
  }
  if (k->sl_size > 0) {
    if (!k->sl_n || !k->sl_k) {
      scs_printf("sum-of-largest-eigenvalues cone arrays missing\n");
      return -1;
    }
    for (i = 0; i < k->sl_size; ++i) {
      if ((k->sl_k[i] >= k->sl_n[i]) || k->sl_k[i] <= 0) {
        scs_printf("sum-of-largest-eigenvalues cone dimension error\n");
        return -1;
      }
    }
  }
#endif
  cone_dims = get_full_cone_dims(k);
  if (cone_dims != d->m) {
    scs_printf("Error: Cone dims %li != rows in A %li\n", (long)cone_dims,
               (long)d->m);
    return -1;
  }
  return 0;
}

/* ================ PSD Block Workspaces and Threads ==================== */

/*
 * PARALLEL PROJECTION OF THE PSD BLOCKS (build option USE_PARALLEL_PSD)
 *
 * What the option does. The real and complex PSD blocks of a cone are
 * projected independently of each other. A build with USE_PARALLEL_PSD, OpenMP
 * and LAPACK shares them out over OpenMP threads, and each thread calls LAPACK
 * for its own blocks. Every other build projects them one after another on
 * the calling thread.
 *
 * Two levels of threads, one at a time. The projection of a block is an
 * eigendecomposition in LAPACK, which runs on BLAS, and a BLAS can have
 * threads of its own. So the cores can be used in two ways:
 * - by the BLAS, inside one block at a time. This is what a build without
 *   the option does, and what a build with it does when it keeps a
 *   projection serial (psd_num_threads): there is no parallel region then,
 *   and the BLAS is as free to use its threads as it is without OpenMP;
 * - by SCS, across the blocks, with the BLAS on one thread inside each.
 * Across the blocks is the faster way when there are several of them: on 4
 * threads, 4 to 16 blocks of n = 256 to 1024 were projected 3 times faster,
 * while the threads of the BLAS inside one block gained nothing below n of
 * about 1000 (OpenBLAS, Accelerate, NVPL). It cannot be faster than its
 * largest block, though, so a cone with one dominant block is left to the
 * BLAS. What must not happen is both at once: t threads of SCS that each
 * start the threads of the BLAS are more threads than cores, and slower than
 * either way alone.
 *
 * What the BLAS does inside the threads of SCS. This decides whether the
 * option pays. SCS cannot ensure it, and cannot even see it for every BLAS,
 * which is why the option is off by default and left to the build, which
 * knows its BLAS:
 * - A BLAS on the OpenMP runtime of SCS (OpenMP builds of OpenBLAS and BLIS,
 *   NVPL, MKL on the matching threading layer) runs on one thread there by
 *   itself. Its parallel regions are nested in ours, and OpenMP gives a
 *   nested region one thread by default. This is the setup the option is
 *   for. The OpenMP thread count (OMP_NUM_THREADS) is then one budget for
 *   both levels: our team during a parallel projection, the threads of the
 *   BLAS everywhere else. It rests on that default: with nested parallelism
 *   switched on (OMP_MAX_ACTIVE_LEVELS above 1), such a BLAS may start
 *   threads inside ours.
 * - A BLAS on another OpenMP runtime (libiomp5 or libomp under a GNU-compiled
 *   SCS) cannot see our region and starts its own threads in each of our
 *   threads. Measured with an OpenBLAS on libomp: 16 threads for 4, at half
 *   the speed. For MKL, each of our threads therefore sets MKL's thread-local
 *   thread count to one for the projection (mkl_set_local_threads). This
 *   goes through a weak reference: it needs ELF and MKL's symbols visible to
 *   SCS.
 * - A BLAS with a thread pool of its own knows nothing of OpenMP. That is a
 *   pthreads OpenBLAS: the default on Debian and Ubuntu, and what Julia and
 *   the SciPy wheels ship. With more than one thread, the BLAS calls of our
 *   threads go to its pool as well and the cores are oversubscribed (measured
 *   up to 3.6x slower than serial), and it has no thread-local setting to
 *   stop that. So when the linked BLAS is a pthreads OpenBLAS on more than
 *   one thread the blocks stay serial (blas_threads_outside_openmp), and
 *   OPENBLAS_NUM_THREADS=1 is what enables the parallel projection. This
 *   check goes through weak references too: it is off for a static OpenBLAS,
 *   one with suffixed symbols, and one loaded privately by another library
 *   (libblastrampoline in Julia). On one thread an OpenBLAS without OpenMP
 *   still makes the BLAS calls of our threads wait for one lock, and the
 *   estimate of the thread count adds that (PSD_PAR_LOCK_COST): small
 *   blocks then stay serial, and the others get fewer threads than there
 *   are. That needs the same symbols.
 * - Accelerate keeps threads of its own as well, which stayed mostly idle
 *   inside ours: at its default settings the parallel projection was 2.3 to
 *   3.3x faster than the serial one on 4 threads, and at most 17% slower
 *   than with Accelerate on one thread (measured with syevr, before the
 *   projection moved to syevd).
 * Neither check exists on Mach-O or Windows, and there is none for other
 * libraries, or for a BLAS other than MKL on another OpenMP runtime. There
 * the BLAS has to be put on one thread for the whole process, or the build
 * left without the option. An environment variable cannot do that for an
 * OpenMP BLAS on another runtime: OMP_NUM_THREADS is read by both runtimes,
 * and an OpenMP build of OpenBLAS ignores OPENBLAS_NUM_THREADS.
 *
 * The split over threads is fixed at setup: blocks in order of decreasing
 * estimated time (psd_block_cost) go to the least loaded thread (LPT
 * scheduling). Thread t projects its blocks (list t) with workspace
 * c->psd_tw[t], whose buffers and LAPACK work arrays are sized for the
 * blocks of list t alone (psd_tw[0] also serves the logdet and
 * sum-of-largest cones, which are projected after the PSD blocks). No
 * workspace is larger than that of the serial projection, so the extra
 * memory is at most (psd_nthreads - 1) times the serial PSD workspace, and
 * less when the large blocks go to few threads. The LAPACK calls for a block
 * pass the work-array sizes of that block (psd_lsz), not the length of the
 * arrays of the workspace: syevd and heevd hand what is left of their work
 * array to sytrd / hetrd and ormtr / unmtr, which lower their block size to
 * fit, so the blocking (and the rounding) can depend on the size passed.
 * Each block then runs the same code with the same sizes on any thread, so
 * with single-threaded BLAS the iterates do not depend on the number of
 * threads.
 *
 * The number of threads is the OpenMP thread count at setup, at most the
 * number of blocks with n > 1. A projection uses fewer threads when the
 * OpenMP thread count has been lowered since: a smaller team still projects
 * every list, each with its own workspace. The number of threads is also
 * bounded by the work. Each thread of a parallel region takes time to wake
 * and to join (PSD_PAR_THREAD_COST), so the region is fastest with fewer
 * threads than there are when the blocks are little work: about
 * sqrt(S / PSD_PAR_THREAD_COST) of them for an estimated serial time S over
 * many blocks. The threads are not used, and the BLAS is left to its own,
 * when
 * - the estimated gain is below PSD_PAR_MIN_GAIN: little work, or one block
 *   with most of it (the parallel loop is no faster than its largest block),
 * - a block larger than PSD_PAR_BLAS_N holds more than half of the work (the
 *   parallel loop then gains less than 2x, while that block alone may run
 *   faster with multithreaded BLAS), or
 * - the workspace is set up inside an active parallel region (the caller's
 *   threads are the parallelism then).
 */
#ifdef PSD_PARALLEL
/* Size of PSD block b: real block b if b < ssize, else complex b - ssize */
static inline scs_int psd_block_n(const ScsCone *k, scs_int b) {
  return b < k->ssize ? k->s[b] : k->cs[b - k->ssize];
}

/* Estimated serial time, in nanoseconds, below which the threads are not
 * used. Zero, because the rules of psd_num_threads keep little work serial
 * already. A large value keeps every projection serial. Those rules take the
 * threads to be awake. Where they sleep between projections (the GNU runtime
 * on Linux, when the rest of an iteration takes long enough), waking them can
 * cost more than a projection of less than a millisecond gains, and a minimum
 * of that size here would avoid it. It is zero there too, on purpose: the
 * same minimum would keep small problems serial, whose threads stay awake and
 * gain the most (an iteration 2.5 times faster, against one up to 27% slower,
 * both measured). OMP_WAIT_POLICY=active, or a build without
 * USE_PARALLEL_PSD, is the switch for the other case. */
#ifndef PSD_PAR_MIN_COST
#define PSD_PAR_MIN_COST (0.)
#endif

/* Time, in nanoseconds, that each thread adds to a parallel region, measured
 * with the GNU runtime on projections that follow each other closely: about
 * 5 us per thread on macOS (Apple M4), where the threads of a team sleep
 * between regions, and 1.5 to 2 us on Linux (Cortex-X925), where they spin
 * for a while after a region. */
#ifndef PSD_PAR_THREAD_COST
#ifdef __APPLE__
#define PSD_PAR_THREAD_COST (5e3)
#else
#define PSD_PAR_THREAD_COST (2e3)
#endif
#endif

/* smallest estimated gain for which the threads are used */
#ifndef PSD_PAR_MIN_GAIN
#define PSD_PAR_MIN_GAIN (1.25)
#endif

/* Largest block for which the threads of the BLAS are not worth more than
 * the loop. With a block of n = 300 or 500 that holds more than half of the
 * work, the loop was 1.6x faster than the serial projection, and the serial
 * projection with a multithreaded BLAS 0.7x to 1.25x (NVPL, OpenBLAS;
 * Accelerate 1.0x at 300). At n = 1000 the two were about equal, and at 2000
 * the BLAS threads were ahead (1.7x to 1.9x against 1.2x to 1.4x). Blocks of
 * n = 1000 to 2000 were 1.3x to 2.1x faster with BLAS threads in #505
 * (OpenBLAS, MKL). All of this was measured with syevr. syevd starts the
 * threads of a multithreaded BLAS at smaller blocks (see EIGENDECOMPOSITION
 * OF THE PSD BLOCKS), which can move the crossover; it has not been measured
 * again. */
#ifndef PSD_PAR_BLAS_N
#define PSD_PAR_BLAS_N (512)
#endif

/* An OpenBLAS that was not built with OpenMP takes one mutex in every BLAS
 * call that needs a work buffer (dsymv: one call per column in the reduction
 * to tridiagonal form), so the calls of the threads of a team wait for each
 * other. A block of size n then takes longer in a team of t threads than
 * alone, by n * (PSD_PAR_LOCK_COST * (t - 1) + PSD_PAR_LOCK_COST2 * (t - 1)^2)
 * nanoseconds. Measured on a Cortex-X925 under Linux (OpenBLAS 0.3.31,
 * pthreads build on one thread), the time per column depended little on the
 * size of the block, for n = 3 to 32: about 0.3, 0.5, 1.1, 1.7 and 4.8 us
 * with 2, 3, 4, 5 and 10 threads. The constants are a third above the fit to
 * that, because the smallest blocks lost more. An OpenMP build did not show
 * this (measured on macOS). These were measured with syevr and not again
 * with syevd; the rule, with the estimate of psd_block_cost, was checked by
 * that estimate's A/B (0.95x the time with this OpenBLAS). */
#ifndef PSD_PAR_LOCK_COST
#define PSD_PAR_LOCK_COST (300.)
#endif
#ifndef PSD_PAR_LOCK_COST2
#define PSD_PAR_LOCK_COST2 (50.)
#endif

/* OpenBLAS: openblas_get_parallel() is 0 sequential, 1 pthreads, 2 OpenMP.
 * MKL: MKL_Set_Num_Threads_Local(n) sets the number of MKL threads for the
 * calling thread only, and returns the previous one (0: the global setting).
 * Weak references, so other BLAS libraries leave them NULL (ELF only: a
 * weak reference to a missing symbol does not link on Mach-O). C linkage,
 * so a C++ build looks up the same symbols. */
#if defined(__GNUC__) && defined(__ELF__)
#ifdef __cplusplus
extern "C" {
#endif
int openblas_get_parallel(void) __attribute__((weak));
int openblas_get_num_threads(void) __attribute__((weak));
int MKL_Set_Num_Threads_Local(int) __attribute__((weak));
#ifdef __cplusplus
}
#endif
#endif

/* Whether BLAS calls would start threads of their own inside OpenMP */
static scs_int blas_threads_outside_openmp(void) {
#if defined(__GNUC__) && defined(__ELF__)
  if (openblas_get_parallel && openblas_get_num_threads) {
    return openblas_get_parallel() == 1 && openblas_get_num_threads() > 1;
  }
#endif
  return 0;
}

/* Whether the BLAS is an OpenBLAS without OpenMP, whose calls from several
 * threads wait for one lock (PSD_PAR_LOCK_COST) */
static scs_int blas_calls_share_lock(void) {
#if defined(__GNUC__) && defined(__ELF__)
  if (openblas_get_parallel) {
    return openblas_get_parallel() != 2;
  }
#endif
  return 0;
}

/* Set the number of MKL threads for the calling thread and return the previous
 * setting, to restore it with. Does nothing, and returns 0, without MKL.
 * SCS_NO_MKL_PSD_GUARD (a test-only build flag) turns this off, so CI can
 * check that the parallel PSD projection oversubscribes without it; the MKL
 * symbol itself cannot be renamed for that, since the MKL linear solver
 * references it too. */
static int mkl_set_local_threads(int nt) {
#if defined(__GNUC__) && defined(__ELF__) && !defined(SCS_NO_MKL_PSD_GUARD)
  if (MKL_Set_Num_Threads_Local) {
    return MKL_Set_Num_Threads_Local(nt);
  }
#endif
  (void)nt;
  return 0;
}

typedef struct {
  scs_float cost;
  scs_int blk;
} PsdBlockCost;

/* Estimated time of the projection of block b, in nanoseconds. It decides
 * the split: blocks go, largest estimate first, to the least loaded thread,
 * and the estimates of all blocks, next to PSD_PAR_THREAD_COST and
 * PSD_PAR_LOCK_COST, decide how many threads are used at all. So only the
 * relative sizes and the order of magnitude matter. n^3 alone would make
 * small blocks far too cheap: a block of n = 10 takes 1/600 to 1/900 of the
 * time of one of n = 300 (the three libraries below), not 1/27000. A block
 * of n = 2 takes LAPACK's short path, and a block of n = 1 is a clamp.
 *
 * Real blocks: a cubic fitted to the times of single blocks with syevd on
 * one thread, as inside the threads of a parallel projection, with NVPL on
 * a Cortex-X925 (DGX Spark, Linux): least squares in log(estimate / time)
 * for n = 3 to 1000. NVPL's times are 0.91-1.17x the estimate, except
 * n = 64 (1.5x). n = 2 is 200 ns, the time with Accelerate and OpenBLAS
 * (NVPL: 1.4 us).
 *
 * The tradeoff. One formula serves every BLAS, and they differ most for
 * small blocks: for n = 3 to 10, Accelerate (Apple M4) and OpenBLAS 0.3.31
 * (pthreads build on one thread, Cortex-X925) take a quarter to two thirds
 * of NVPL's time, and Accelerate grows more slowly for large blocks (55 ms
 * at n = 1000, NVPL 116 ms). Too low an estimate keeps blocks serial that
 * threads would project faster; too high gives threads to projections that
 * are faster on one. Following NVPL, the highest times, was faster with two
 * of the three libraries. On the block structures whose split changed
 * (#507's 82 and the SDPs of SDPLIB and the benchmark suite with several
 * blocks, 2 to 4 threads), the projection took, against the previous
 * estimate (below), 0.86x the time with NVPL (4 blocks of n = 5 and 6 of
 * n = 3 now take threads, 1.6-2.8x faster), 0.95x with OpenBLAS (lock
 * counted, PSD_PAR_LOCK_COST) and 1.05x with Accelerate, where projections
 * of 9 to 16 us made of 8 to 16 blocks of n = 3 to 6 now take threads and
 * are up to 2.3x slower. A fit through NVPL and Accelerate together
 * (0.5-0.8x this estimate for n = 3 to 10) took 1.10x with NVPL, 1.18x with
 * OpenBLAS and 1.00x with Accelerate. NVPL's 1.4 us for n = 2 took 0.4-0.7x
 * with NVPL on 16 to 64 such blocks, but up to 1.8x with the other two.
 *
 * Complex blocks: the estimate from before syevd (#543), fitted to heevr's
 * times with Accelerate and OpenBLAS. heevd with NVPL takes 0.9-1.4x of it
 * for n = 3 to 100 and 1.2-1.55x from n = 150 (Accelerate 1.2-1.3x,
 * OpenBLAS up to 1.85x at n = 500). The cubic fitted to heevd with NVPL
 * (1170 + 345 n + 60 n^2 + 0.45 n^3) was slower: 1.20x with NVPL, 1.09x with
 * OpenBLAS and 1.08x with Accelerate on the complex structures whose split
 * changed, mostly because one block of n = 300 with four of n = 100 then
 * looked like a gain of 1.24x, below PSD_PAR_MIN_GAIN, and stayed serial
 * (1.2-1.55x slower).
 *
 * Previous estimate, fitted to syevr: 600 + 250 n + 45 n^2 + 0.06 n^3. */
static scs_float psd_block_cost(const ScsCone *k, scs_int b) {
  scs_float n = (scs_float)psd_block_n(k, b);
  if (b < k->ssize)
    return n < 3. ? (n < 2. ? 2. : 200.)
                  : 940. + n * (310. + n * (45. + 0.055 * n));
  return n < 3. ? (n < 2. ? 2. : 300.)
                : 600. + n * (250. + n * (70. + 0.24 * n));
}

/* decreasing cost, ties by block index so the split is deterministic */
static int cmp_psd_block_cost(const void *a, const void *b) {
  const PsdBlockCost *x = (const PsdBlockCost *)a;
  const PsdBlockCost *y = (const PsdBlockCost *)b;
  if (x->cost != y->cost)
    return x->cost < y->cost ? 1 : -1;
  return (x->blk > y->blk) - (x->blk < y->blk);
}

/* LPT: each block, largest first, to the least loaded thread; each thread's
 * list largest first. Fills psd_tptr and psd_blk for psd_nthreads > 1. */
static scs_int assign_psd_blocks(ScsConeWork *c, const ScsCone *k) {
  scs_int nblk = k->ssize + k->cssize, nt = c->psd_nthreads, i, j, t;
  scs_int status = 0, *pos = (scs_int *)scs_calloc(nblk, sizeof(scs_int));
  scs_float *load = (scs_float *)scs_calloc(nt, sizeof(scs_float));
  PsdBlockCost *bc = (PsdBlockCost *)scs_calloc(nblk, sizeof(PsdBlockCost));
  if (!pos || !load || !bc) {
    status = -1;
    goto out;
  }
  for (i = 0; i < nblk; ++i) {
    bc[i].blk = i;
    bc[i].cost = psd_block_cost(k, i);
  }
  qsort(bc, nblk, sizeof(PsdBlockCost), cmp_psd_block_cost);
  memset(c->psd_tptr, 0, (nt + 1) * sizeof(scs_int));
  for (j = 0; j < nblk; ++j) {
    scs_int tmin = 0;
    for (t = 1; t < nt; ++t) {
      if (load[t] < load[tmin])
        tmin = t;
    }
    load[tmin] += bc[j].cost;
    pos[j] = tmin;
    c->psd_tptr[tmin + 1]++;
  }
  for (t = 0; t < nt; ++t)
    c->psd_tptr[t + 1] += c->psd_tptr[t];
  /* psd_tptr[t] is the fill cursor of thread t, then shift back to list
   * starts */
  for (j = 0; j < nblk; ++j)
    c->psd_blk[c->psd_tptr[pos[j]]++] = bc[j].blk;
  for (t = nt; t > 0; --t)
    c->psd_tptr[t] = c->psd_tptr[t - 1];
  c->psd_tptr[0] = 0;
out:
  if (pos)
    scs_free(pos);
  if (load)
    scs_free(load);
  if (bc)
    scs_free(bc);
  return status;
}

/* Estimated time of the parallel loop with nt threads over nblk blocks, which
 * take total in all and cmax for the largest. With a lock in the BLAS
 * (PSD_PAR_LOCK_COST), ncol is the number of columns of these blocks and nmax
 * that of the largest; without one, both are zero. */
static scs_float psd_par_time(scs_int nt, scs_int nblk, scs_float total,
                              scs_float cmax, scs_float ncol, scs_float nmax) {
  scs_float lock =
      (nt - 1) * (PSD_PAR_LOCK_COST + (nt - 1) * PSD_PAR_LOCK_COST2);
  scs_float list =
      (scs_float)((nblk + nt - 1) / nt) * ((total + lock * ncol) / nblk);
  return MAX(cmax + lock * nmax, list) + nt * PSD_PAR_THREAD_COST;
}
#endif

/* Number of threads for the PSD blocks (see above) */
static scs_int psd_num_threads(const ScsCone *k) {
#ifdef PSD_PARALLEL
  scs_int b, t, ntmax, nbig = 0, nmax = 0, ncol = 0;
  scs_int nt = omp_get_max_threads(), lock = blas_calls_share_lock();
  scs_float total = 0., cmax = 0., est, best = 0., lcol, lmax;
  for (b = 0; b < k->ssize + k->cssize; ++b) {
    scs_int n = psd_block_n(k, b);
    scs_float cost = psd_block_cost(k, b);
    total += cost;
    cmax = MAX(cmax, cost);
    nmax = MAX(nmax, n);
    if (n > 1) {
      nbig++;
      ncol += n;
    }
  }
  nt = MIN(nt, nbig);
  if (nt < 2 || total < PSD_PAR_MIN_COST || omp_in_parallel() ||
      blas_threads_outside_openmp())
    return 1;
  if (2 * cmax > total && nmax > PSD_PAR_BLAS_N)
    return 1;
  /* With t threads the parallel loop takes about as long as its longest
   * list, which is at least the largest block and at least its share of the
   * blocks at their average time, plus the time of the t threads. Find the
   * fewest threads with the smallest such estimate. Without a lock in the
   * BLAS, then take the most threads that come within 5% of it: the estimate
   * is flat where one block bounds the loop, and more threads are then safer
   * for the lists of the other blocks. With a lock every thread slows the
   * others down. */
  lcol = lock ? (scs_float)ncol : 0.;
  lmax = lock ? (scs_float)nmax : 0.;
  ntmax = nt;
  for (t = 2; t <= ntmax; ++t) {
    est = psd_par_time(t, nbig, total, cmax, lcol, lmax);
    if (t == 2 || est < best) {
      best = est;
      nt = t;
    }
  }
  for (t = ntmax; !lock && t > nt; --t) {
    if (psd_par_time(t, nbig, total, cmax, lcol, lmax) <= 1.05 * best)
      nt = t; /* the last pass: t > nt no longer holds */
  }
  return total < PSD_PAR_MIN_GAIN * best ? 1 : nt;
#else
  (void)k;
  return 1;
#endif
}

/* Offsets of the PSD blocks and their split over threads */
static scs_int set_up_psd_schedule(ScsConeWork *c, const ScsCone *k) {
  scs_int nblk = k->ssize + k->cssize, nt = psd_num_threads(k), b;
  c->psd_tptr = (scs_int *)scs_calloc(nt + 1, sizeof(scs_int));
  c->psd_blk = (scs_int *)scs_calloc(nblk, sizeof(scs_int));
  c->psd_off = (scs_int *)scs_calloc(nblk, sizeof(scs_int));
  if (!c->psd_tptr || !c->psd_blk || !c->psd_off)
    return -1;
  c->psd_nthreads = nt;
  c->psd_len = 0;
  for (b = 0; b < nblk; ++b) {
    c->psd_off[b] = c->psd_len;
    c->psd_len += b < k->ssize ? get_sd_cone_size(k->s[b])
                               : get_csd_cone_size(k->cs[b - k->ssize]);
    c->psd_blk[b] = b;
  }
  c->psd_tptr[1] = nblk;
#ifdef PSD_PARALLEL
  if (nt > 1)
    return assign_psd_blocks(c, k);
#endif
  return 0;
}

#ifdef USE_LAPACK
/* Buffers of a workspace for eigenvalues up to n, real matrices up to nr,
 * complex ones up to nc, and results of the spectral cones up to nz (none
 * if 0). No buffer for the eigenvectors: syevd and heevd return them in Xs
 * and cXs, and the PSD projections put their results in the work arrays. */
static scs_int alloc_psd_bufs(ScsPsdWork *w, blas_int n, blas_int nr,
                              blas_int nc, blas_int nz) {
  w->e = (scs_float *)scs_calloc(n, sizeof(scs_float));
  w->gamma = (scs_float *)scs_calloc(n, sizeof(scs_float));
  if (!w->e || !w->gamma)
    return -1;
  if (nr > 0 &&
      !(w->Xs = (scs_float *)scs_calloc((size_t)nr * nr, sizeof(scs_float))))
    return -1;
  if (nz > 0 &&
      !(w->Z = (scs_float *)scs_calloc((size_t)nz * nz, sizeof(scs_float))))
    return -1;
  if (nc > 0 && !(w->cXs = (scs_complex_float *)scs_calloc(
                      (size_t)nc * nc, sizeof(scs_complex_float))))
    return -1;
  return 0;
}

/* LAPACK work arrays of a workspace, of the given lengths */
static scs_int alloc_psd_work_arrays(ScsPsdWork *w, blas_int lwork,
                                     blas_int liwork, blas_int lcwork) {
  w->lwork = lwork;
  w->liwork = liwork;
  w->lcwork = lcwork;
  if (lwork > 0 &&
      !(w->work = (scs_float *)scs_calloc(lwork, sizeof(scs_float))))
    return -1;
  if (liwork > 0 &&
      !(w->iwork = (blas_int *)scs_calloc(liwork, sizeof(blas_int))))
    return -1;
  if (lcwork > 0 && !(w->cwork = (scs_complex_float *)scs_calloc(
                          lcwork, sizeof(scs_complex_float))))
    return -1;
  return 0;
}

/*
 * EIGENDECOMPOSITION OF THE PSD BLOCKS
 *
 * The projection of a PSD block is an eigendecomposition by LAPACK's syevd
 * (heevd for a complex block): divide and conquer, with the eigenvectors
 * returned in place of the block. SCS used syevr / heevr (MRRR) before.
 * - Speed, real blocks. On SDPLIB, syevd took 31-51% less cone time per
 *   iteration than syevr, and the iterations 28-47% less time in all
 *   (geometric means: NVPL and OpenBLAS on a 20-core Arm GB10 at 1 and 4
 *   threads on its fast cores and at 20 threads, Accelerate on an Apple M4);
 *   with MKL on an AMD Zen 3 at 1 and 4 threads, 29-38% less cone time on 14
 *   of the problems. The gain is largest on blocks with exactly repeated
 *   eigenvalues, which MRRR is slow on: 2.7-2.9x on thetaG11. With the BLAS
 *   on 1 or 4 fast cores of the GB10 no problem was more than 2.3% slower
 *   (qap8, NVPL; OpenBLAS on one thread on full solves: a run of qpG11 capped
 *   at 75 iterations was 19% slower), nor more than 1.8% with Accelerate (82
 *   problems), and none was slower with MKL. On one of the GB10's slow cores
 *   (Cortex-A725) the gain was about 20% on 10 problems with NVPL, and qap8
 *   was 11-16% slower (two checks of 3 runs on a loaded machine).
 * - Speed, complex blocks (none in SDPLIB). heevd rather than heevr (with
 *   the heev fallback), because heevr did not save the 3% of solve time
 *   that was set, before measuring, as the bar for a second routine. On 28
 *   problems with complex blocks of n = 2 to 300 (quantum fidelity, PPT,
 *   channel capacity, PhaseLift, phase synchronization, complex max-cut and
 *   the test suite's) at eps 1e-4 and 1e-6, with Accelerate on an Apple M4
 *   (1 and default threads), MKL on an AMD Zen 3 (1 and 4 threads), and
 *   NVPL and OpenBLAS on one fast core of the GB10, both solved all 324
 *   runs, in the same iterations, and heevr's shifted geometric mean of the
 *   solve time (shift 1 s) was 0.997 of heevd's. heevr took more cone time
 *   per iteration with every LAPACK: 1.07x with Accelerate, 1.15-1.21x with
 *   MKL, 1.16x with NVPL and 1.25x with OpenBLAS; with MKL more at every
 *   block size (1.10-1.16x above n = 150), and with every LAPACK more for
 *   n <= 64 (1.07-1.43x). It took less only on larger blocks, with
 *   Accelerate (0.91x for n = 65 to 150 and 0.71-0.74x above, so that the
 *   solves with blocks of n = 200 to 300 took 0.67-0.83x as long) and with
 *   NVPL above n = 150 (0.91x). Where about half of the eigenvalues came in
 *   pairs closer than 1e-8 relative to each other (the PPT problems), heevr
 *   took 1.07-1.69x as long as heevd to solve, with every LAPACK.
 * - Accuracy. syevd was as accurate or more (errors of 2-4e-15 relative to the
 *   largest entry of the projection against 3-9e-15, on blocks with exact
 *   multiplets).
 * - Memory. syevd needs a work array of 1 + 6n + 2n^2 floats (heevd: n^2 + 2n
 *   complex numbers and 1 + 5n + 2n^2 floats), where syevr needed about 26n
 *   and a separate n x n array for the eigenvectors, which syevd does not.
 *   The work array also takes the result of the reconstruction (syrk / herk,
 *   n^2 entries), so the workspace of a real block grows from about 2n^2 to
 *   3n^2 floats, that of a complex block from 4n^2 to 6n^2: +n^2 floats per
 *   real block, +8 MB at n = 1000, +200 MB at n = 5000 (+266 MB measured in
 *   peak memory) and +390 MB at n = 7000 (+532 MB measured), in double.
 * - Threads. syevd starts the threads of a multithreaded BLAS at smaller
 *   blocks than syevr. With the BLAS on all 20 cores of the GB10 (10 fast and
 *   10 slow), blocks of n = 45 to 130 took up to 4x longer, so that control7
 *   took 2.1x the cone time per iteration (NVPL) and qap7 1.6x (OpenBLAS), and
 *   OpenBLAS kept up to 20 cores busy where syevr kept 1.4 to 12. Small
 *   blocks are fastest with the BLAS on one thread, with either routine.
 * The fallback (#531, #532): when syevd / heevd reports no convergence (info
 * > 0), the block is decomposed again by syev / heev (QR iteration: 1.4-3.6x
 * the cone time of syevd when forced on 11 SDPLIB problems, which it solved
 * as the other builds did). Apple Accelerate's dsyevd (classic LAPACK
 * interface) failed with info = 322202 on a thetaG11 iterate (n = 801) that
 * syev decomposes, and its dsyevr on the a*I + b*11' block of #531. The
 * test-only build flag SCS_FORCE_PSD_FALLBACK sends every block through the
 * fallback, so CI covers it.
 *
 * Work-array sizes: LAPACK's workspace query at the size of each block
 * (psd_lsz), and at least the minimum that LAPACK documents (psd_min_sizes).
 * Every call for a block passes these sizes, see PARALLEL PROJECTION above.
 * With 32-bit LAPACK integers the minimum fits up to n = 32766: setup fails
 * with a message for a larger block, before any LAPACK call (LAPACK computes
 * the sizes in its own integers, which would overflow).
 */

/* The smallest work arrays that syevd and heevd accept for a block of size
 * n >= 2 with eigenvectors (jobz = "V"), from the LAPACK documentation:
 *   syevd: lwork >= 1 + 6n + 2n^2, liwork >= 3 + 5n,
 *   heevd: lrwork >= 1 + 5n + 2n^2, lwork >= 2n + n^2 (complex numbers),
 *          liwork >= 3 + 5n.
 * In lw (work of a real block, rwork of a complex one), lcw (complex work, 0
 * for a real block) and liw. In double, which holds them exactly far beyond
 * where they overflow a blas_int. lw is the largest. */
static void psd_min_sizes(double n, scs_int cplx, double *lw, double *lcw,
                          double *liw) {
  *lw = cplx ? 1. + 5. * n + 2. * n * n : 1. + 6. * n + 2. * n * n;
  *lcw = cplx ? 2. * n + n * n : 0.;
  *liw = 3. + 5. * n;
}

/* Whether the work arrays of a PSD block of size n >= 2 fit in a blas_int.
 * Prints why not if they do not. */
static scs_int psd_block_fits(scs_int n, scs_int cplx) {
  double lw, lcw, liw;
  psd_min_sizes((double)n, cplx, &lw, &lcw, &liw);
  if (lw <= (double)SCS_BLAS_INT_MAX)
    return 1;
  scs_printf("FATAL: a %sPSD cone of size %li needs a LAPACK work array of "
             "%.0f entries, more than LAPACK integers can count here (%.0f). "
             "Use SCS with a BLAS and LAPACK with 64-bit integers (build "
             "option BLAS64) for cones this large, and with 64-bit SCS "
             "integers (DLONG) as well above size 46340.\n",
             cplx ? "complex " : "", (long)n, lw, (double)SCS_BLAS_INT_MAX);
  return 0;
}

/* A size that LAPACK returns as a floating-point number, rounded up by one
 * (as SCS always has: a float cannot hold every size above 2^24), or the
 * minimum min_size when that is smaller or does not fit in a blas_int */
static blas_int lapack_size(double queried, double min_size) {
  double sz = queried + 1.;
  /* also catches NaN */
  if (!(sz >= min_size && sz <= (double)SCS_BLAS_INT_MAX))
    sz = min_size;
  return (blas_int)sz;
}

/* Work-array sizes of a PSD block of size n >= 2 (complex if cplx), in sz,
 * from LAPACK's workspace query. w supplies the matrix and eigenvalue
 * arguments, which a query does not touch: w->Xs (w->cXs for a complex
 * block) and w->e must hold a block of size n. */
static scs_int query_psd_block_sizes(ScsPsdWork *w, blas_int n, scs_int cplx,
                                     ScsPsdLapackSizes *sz) {
  blas_int neg_one = -1, info = 0, liwork_opt = 0;
  double lw_min, lcw_min, liw_min;
  psd_min_sizes((double)n, cplx, &lw_min, &lcw_min, &liw_min);
  if (cplx) {
    scs_complex_float lcwork_opt = {0.0};
    scs_float lrwork_opt = 0.0;
    BLASC(heevd)("V", "L", &n, SCS_BLAS_COMPLEX_CAST(w->cXs), &n, w->e,
                 SCS_BLAS_COMPLEX_CAST(&lcwork_opt), &neg_one, &lrwork_opt,
                 &neg_one, &liwork_opt, &neg_one, &info);
    sz->lwork = lapack_size((double)lrwork_opt, lw_min);
    sz->lcwork = lapack_size((double)lcwork_opt[0], lcw_min);
  } else {
    scs_float lwork_opt = 0.0;
    BLAS(syevd)("V", "L", &n, w->Xs, &n, w->e, &lwork_opt, &neg_one,
                &liwork_opt, &neg_one, &info);
    sz->lwork = lapack_size((double)lwork_opt, lw_min);
    sz->lcwork = 0;
  }
  if (info != 0) {
    scs_printf("FATAL: %s workspace query failure, info = %li\n",
               cplx ? "heevd" : "syevd", (long)info);
    return -1;
  }
  sz->liwork = MAX(liwork_opt, (blas_int)liw_min);
  return 0;
}

/* Sets up workspace t (c->psd_tw[t]) for the PSD blocks of list t: buffers
 * for its largest real and complex blocks, the LAPACK work-array sizes of
 * each of its blocks (c->psd_lsz) and work arrays for the largest of them.
 * The workspace queries pass this workspace's buffers as the matrix and
 * eigenvalue arguments, which a query does not touch; consecutive blocks of
 * the same size and kind share one query, and blocks of size 1 call no
 * LAPACK and keep zero sizes. n_spec is the largest logdet or
 * sum-of-largest-eigenvalues cone that this workspace also serves (only
 * workspace 0 does), or 0: it then gets the buffers of that size and the
 * work array of the size those projections pass to syev (c->lwork_syev). */
static scs_int set_up_psd_work(ScsConeWork *c, const ScsCone *k, scs_int t,
                               blas_int n_spec) {
  ScsPsdWork *w = &(c->psd_tw[t]);
  /* the last query, real and complex */
  ScsPsdLapackSizes last[2] = {{0, 0, 0}, {0, 0, 0}};
  blas_int last_n[2] = {0, 0};
  blas_int nr = 0, nc = 0, lwork = 0, liwork = 0, lcwork = 0;
  scs_int j, j0 = 0, j1 = 0; /* list t: psd_blk[j0 .. j1 - 1] */
  if (c->psd_nthreads > 0) {
    j0 = c->psd_tptr[t];
    j1 = c->psd_tptr[t + 1];
  }
  for (j = j0; j < j1; ++j) {
    scs_int b = c->psd_blk[j];
    if (b < k->ssize)
      nr = MAX(nr, (blas_int)k->s[b]);
    else
      nc = MAX(nc, (blas_int)k->cs[b - k->ssize]);
  }
  /* a block of size 1 is a clamp and needs no buffer */
  nr = nr > 1 ? nr : 0;
  nc = nc > 1 ? nc : 0;
  if (alloc_psd_bufs(w, MAX(MAX(MAX(nr, nc), n_spec), 1), MAX(nr, n_spec), nc,
                     n_spec) < 0)
    return -1;
  for (j = j0; j < j1; ++j) {
    scs_int b = c->psd_blk[j], cplx = b >= k->ssize;
    blas_int n = (blas_int)(cplx ? k->cs[b - k->ssize] : k->s[b]);
    if (n < 2)
      continue;
    if (n != last_n[cplx]) {
      if (query_psd_block_sizes(w, n, cplx, &last[cplx]) < 0)
        return -1;
      last_n[cplx] = n;
    }
    c->psd_lsz[b] = last[cplx];
    lwork = MAX(lwork, last[cplx].lwork);
    liwork = MAX(liwork, last[cplx].liwork);
    lcwork = MAX(lcwork, last[cplx].lcwork);
  }
#ifdef USE_SPECTRAL_CONES
  /* The logdet and sum-of-largest projections eigendecompose with syev in
   * this workspace's Xs, e and work. They pass c->lwork_syev, syev's own
   * workspace query for the largest of these cones, not the length of the
   * work array: that length depends on the PSD blocks of list 0, and so on
   * the number of threads, which the blocking inside syev (and its rounding)
   * then would too. Neither size bounds the other: a spectral cone can be
   * larger than every PSD block, and syev's optimum (34n with Accelerate and
   * netlib) is above syevd's minimum for n < 14. */
  if (n_spec > 0) {
    blas_int neg_one = -1, info = 0;
    scs_float wkopt = 0.0;
    BLAS(syev)("V", "L", &n_spec, w->Xs, &n_spec, w->e, &wkopt, &neg_one,
               &info);
    if (info != 0) {
      scs_printf("FATAL: syev workspace query failure, info = %li\n",
                 (long)info);
      return -1;
    }
    /* syev's documented minimum: max(1, 3n - 1) */
    c->lwork_syev =
        lapack_size((double)wkopt, MAX(1., 3. * (double)n_spec - 1.));
    lwork = MAX(lwork, c->lwork_syev);
  }
#endif
  return alloc_psd_work_arrays(w, lwork, liwork, lcwork);
}
#endif

/* ======================== Workspace Setup ============================= */

/*
 * Consolidated setup for Real PSD, Complex PSD, and Spectral cones.
 */
static scs_int set_up_cone_work_spaces(ScsConeWork *c, const ScsCone *k) {
  scs_int i;
#ifdef USE_LAPACK
  scs_int t;
  /* the largest logdet or sum-of-largest-eigenvalues cone, which workspace 0
   * serves too (0 without one) */
  blas_int n_spec = 0;

  /* Fail before any allocation or LAPACK call if a PSD block's work arrays
   * do not fit in a blas_int */
  for (i = 0; i < k->ssize; ++i) {
    if (k->s[i] > 1 && !psd_block_fits(k->s[i], 0))
      return -1;
  }
  for (i = 0; i < k->cssize; ++i) {
    if (k->cs[i] > 1 && !psd_block_fits(k->cs[i], 1))
      return -1;
  }

#ifdef USE_SPECTRAL_CONES
  /* LAPACK Query Variables */
  blas_int neg_one = -1, info = 0;
  scs_float wkopt = 0.0;

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

  /* The logdet and sum-of-largest projections work in workspace 0 */
  if (k->dsize > 0 || k->sl_size > 0)
    n_spec = MAX(n_max_logdet, n_max_sl);
#endif

  /* 1-2. One eigendecomposition workspace per PSD thread, each set up for
   * the blocks of its own list with syevd's and heevd's workspace queries
   * (see EIGENDECOMPOSITION OF THE PSD BLOCKS), and workspace 0 for the
   * logdet and sum-of-largest cones as well. heevd's real 'rwork' array is
   * the real 'work' array of a workspace. */
  c->psd_tw = (ScsPsdWork *)scs_calloc(MAX(c->psd_nthreads, 1),
                                       sizeof(ScsPsdWork));
  c->psd_lsz = (ScsPsdLapackSizes *)scs_calloc(MAX(k->ssize + k->cssize, 1),
                                               sizeof(ScsPsdLapackSizes));
  if (!c->psd_tw || !c->psd_lsz)
    return -1;
  for (t = 0; t < MAX(c->psd_nthreads, 1); ++t) {
    if (set_up_psd_work(c, k, t, t == 0 ? n_spec : 0) < 0)
      return -1;
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

/* ====================== Cone Projections ============================= */

/*
 * Projection: Real Semi-Definite Cone
 */
#ifdef USE_LAPACK
/* Unpack the lower triangle of an n x n block, stored column-wise in X, into
 * the full column-major buffer Xs (only the lower triangle is written). */
static void unpack_sd_lower(const scs_float *X, scs_int n, scs_float *Xs) {
  scs_int i;
  for (i = 0; i < n; ++i) {
    memcpy(&(Xs[i * (n + 1)]), &(X[i * n - ((i - 1) * i) / 2]),
           (n - i) * sizeof(scs_float));
  }
}

/* Unpack a complex n x n block (lower triangle, packed as reals in X) into
 * the full column-major complex buffer cXs. */
static void unpack_csd_lower(const scs_float *X, scs_int n,
                             scs_complex_float *cXs) {
  scs_int i;
  for (i = 0; i < n - 1; ++i) {
    cXs[i * (n + 1)][0] = X[i * (2 * n - i)]; /* Diagonal (Real) */
    cXs[i * (n + 1)][1] = 0.0;
    memcpy(&(cXs[i * (n + 1) + 1]), &(X[i * (2 * n - i) + 1]),
           2 * (n - i - 1) * sizeof(scs_float));
  }
  cXs[n * n - 1][0] = X[n * n - 1]; /* Last element */
  cXs[n * n - 1][1] = 0.0;
}

#ifdef SCS_FORCE_PSD_FALLBACK
/* Test-only build flag: every syevd / heevd call of the PSD projection is
 * treated as if it had not converged, after it has run and overwritten the
 * block, so a test build projects every PSD block through the syev / heev
 * fallback, which no LAPACK we test with reaches on the tests otherwise.
 * psd_fallbacks counts the fallbacks, so that the tests can check that the
 * flag reached the compiler. */
scs_int SCS(psd_fallbacks) = 0;
#endif

/* Whether to decompose a block again with syev / heev after syevd / heevd
 * returned info: when it reports no convergence (info > 0). Always, unless
 * info < 0, under SCS_FORCE_PSD_FALLBACK. */
static scs_int psd_fall_back(blas_int info) {
#ifdef SCS_FORCE_PSD_FALLBACK
  if (info < 0)
    return 0;
#ifdef _OPENMP
#pragma omp atomic
#endif
  SCS(psd_fallbacks)++;
  return 1;
#else
  return info > 0;
#endif
}
#endif

/* Project one svec-packed block onto the PSD cone. When rho (the
 * block's slice of the diagonal metric R) is non-uniform it must have
 * the rank-one structure w_ij = delta_i * delta_j (enforced by
 * fit_psd_block_weights). The projection under the induced R^{-1}
 * metric is then a diagonal congruence sandwich: with gamma_i =
 * w_ii^{-1/4}, project Gamma X Gamma and unscale. Positive homogeneity
 * of the cone makes any uniform factor in rho drop out, so the
 * historical scalar-per-block case reduces to the unscaled path. sz holds
 * the LAPACK work-array sizes of the block (psd_lsz), which w has room for. */
static scs_int proj_semi_definite_cone(scs_float *X, const scs_int n,
                                       ScsPsdWork *w,
                                       const ScsPsdLapackSizes *sz,
                                       const scs_float *rho) {
  scs_int scaled = 0;
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
  blas_int lwork = sz->lwork, liwork = sz->liwork;
  scs_float zero = 0., one = 1., sqrt2 = SQRTF(2.0), sqrt2_inv = 1.0 / sqrt2;
  scs_float sq_eig;
  scs_int first_idx = -1;

  /* Diagonal-congruence pre-scale (skip when the metric is uniform on
   * this block so the historical path stays bit-identical). */
  if (rho) {
    scs_float wmin = rho[0], wmax = rho[0];
    for (i = 1; i < n; ++i) {
      scs_float wii = rho[i * n - ((i - 1) * i) / 2];
      wmin = MIN(wmin, wii);
      wmax = MAX(wmax, wii);
    }
    if (wmax > wmin * (1. + 1e-12)) {
      scs_int jj, ii;
      scaled = 1;
      for (i = 0; i < n; ++i) {
        w->gamma[i] = POWF(rho[i * n - ((i - 1) * i) / 2], -0.25);
      }
      for (jj = 0; jj < n; ++jj) {
        scs_float gj = w->gamma[jj];
        scs_float *col = &(X[jj * n - ((jj - 1) * jj) / 2]);
        for (ii = jj; ii < n; ++ii) {
          col[ii - jj] *= w->gamma[ii] * gj;
        }
      }
    }
  }

  /* Copy lower triangular part to full matrix buffer Xs, scale diagonals by
   * sqrt(2) */
  unpack_sd_lower(X, n, w->Xs);
  BLAS(scal)(&nb, &sqrt2, w->Xs, &nb_p1);

  /* Solve Eigenproblem: Xs = Z * diag(e) * Z', with syevd, which returns the
   * eigenvectors Z in Xs (see EIGENDECOMPOSITION OF THE PSD BLOCKS) */
  BLAS(syevd)("V", "L", &nb, w->Xs, &nb, w->e, w->work, &lwork, w->iwork,
              &liwork, &info);
  if (psd_fall_back(info)) {
    /* syevd reported no convergence, as Apple Accelerate's did on a thetaG11
     * iterate. syevd has overwritten Xs, so unpack again and retry with the
     * QR-based syev, which also returns the eigenvectors in Xs. syev needs
     * lwork >= 3n - 1, well below syevd's minimum of 1 + 6n + 2n^2, so the
     * work array is large enough. */
    unpack_sd_lower(X, n, w->Xs);
    BLAS(scal)(&nb, &sqrt2, w->Xs, &nb_p1);
    BLAS(syev)("V", "L", &nb, w->Xs, &nb, w->e, w->work, &lwork, &info);
  }
  /* info < 0: bad argument. info > 0: no convergence (from the fallback too).
   * Either way there is no eigendecomposition to project with, so fail (the
   * caller reports it: this may run in a PSD thread, which must not print). */
  if (info != 0)
    return -1;

  /* Filter negative eigenvalues and scale eigenvectors */
  /* Note: e is in ascending order */
  for (i = 0; i < n; ++i) {
    if (w->e[i] > 0) {
      if (first_idx == -1) {
        first_idx = i;
      }
      sq_eig = SQRTF(w->e[i]);
      BLAS(scal)(&nb, &sq_eig, &w->Xs[i * n], &one_int);
    }
  }
  if (first_idx == -1) {
    /* there are no positive eigenvalues, set X to 0 and return */
    memset(X, 0, sizeof(scs_float) * get_sd_cone_size(n));
    return 0;
  }
  /* Reconstruct Z * Z' (Z, scaled, in Xs) in the work array, which has room
   * for it: lwork >= 1 + 6n + 2n^2 (psd_min_sizes) */
  ncols_z = (blas_int)(n - first_idx);
  BLAS(syrk)("Lower", "NoTrans", &nb, &ncols_z, &one, &w->Xs[first_idx * n],
             &nb, &zero, w->work, &nb);

  /* Rescale diagonals by 1/sqrt(2) */
  BLAS(scal)(&nb, &sqrt2_inv, w->work, &nb_p1);

  /* Extract lower triangular matrix back to X */
  for (i = 0; i < n; ++i) {
    memcpy(&(X[i * n - ((i - 1) * i) / 2]), &(w->work[i * (n + 1)]),
           (n - i) * sizeof(scs_float));
  }

  /* Undo the congruence pre-scale */
  if (scaled) {
    scs_int jj, ii;
    for (jj = 0; jj < n; ++jj) {
      scs_float gj = w->gamma[jj];
      scs_float *col = &(X[jj * n - ((jj - 1) * jj) / 2]);
      for (ii = jj; ii < n; ++ii) {
        col[ii - jj] /= w->gamma[ii] * gj;
      }
    }
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
                                               ScsPsdWork *w,
                                               const ScsPsdLapackSizes *sz) {
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
  blas_int lwork = sz->lwork, lcwork = sz->lcwork, liwork = sz->liwork;
  scs_float zero = 0., one = 1.;
  scs_int first_idx = -1;

  /* Complex constants */
  scs_complex_float csqrt2 = {0.0}, csqrt2_inv = {0.0}, csq_eig = {0.0};
  csqrt2[0] = SQRTF(2.0);
  csqrt2_inv[0] = 1.0 / csqrt2[0];

  /* Unpack X (real array) into cXs (complex matrix), scale diagonals by
   * sqrt(2) */
  unpack_csd_lower(X, n, w->cXs);
  BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2),
              SCS_BLAS_COMPLEX_CAST(w->cXs), &nb_p1);

  /* Solve Eigenproblem with heevd, which returns the eigenvectors in cXs.
   * Note: w->work acts as rwork here. */
  BLASC(heevd)("V", "L", &nb, SCS_BLAS_COMPLEX_CAST(w->cXs), &nb, w->e,
               SCS_BLAS_COMPLEX_CAST(w->cwork), &lcwork, w->work, &lwork,
               w->iwork, &liwork, &info);
  if (psd_fall_back(info)) {
    /* heevd reported no convergence: unpack again and retry with heev, as
     * the real case retries with syev. heev needs lwork >= 2n - 1 complex
     * numbers and rwork >= 3n - 2, below heevd's minima of 2n + n^2 and
     * 1 + 5n + 2n^2. */
    unpack_csd_lower(X, n, w->cXs);
    BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2),
                SCS_BLAS_COMPLEX_CAST(w->cXs), &nb_p1);
    BLASC(heev)("V", "L", &nb, SCS_BLAS_COMPLEX_CAST(w->cXs), &nb, w->e,
                SCS_BLAS_COMPLEX_CAST(w->cwork), &lcwork, w->work, &info);
  }
  /* info != 0: no eigendecomposition, so fail (see the real case) */
  if (info != 0)
    return -1;

  /* Reconstruct */
  for (i = 0; i < n; ++i) {
    if (w->e[i] > 0) {
      if (first_idx == -1) {
        first_idx = i;
      }
      csq_eig[0] = SQRTF(w->e[i]);
      BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csq_eig),
                  SCS_BLAS_COMPLEX_CAST(&w->cXs[i * n]), &one_int);
    }
  }
  if (first_idx == -1) {
    /* there are no positive eigenvalues, set X to 0 and return */
    memset(X, 0, sizeof(scs_float) * get_csd_cone_size(n));
    return 0;
  }

  /* cwork = Z * Z' (Z, scaled, in cXs): the complex work array has room for
   * it, lcwork >= 2n + n^2 (psd_min_sizes) */
  ncols_z = (blas_int)(n - first_idx);
  BLASC(herk)("Lower", "NoTrans", &nb, &ncols_z, &one,
              SCS_BLAS_COMPLEX_CAST(&w->cXs[first_idx * n]), &nb, &zero,
              SCS_BLAS_COMPLEX_CAST(w->cwork), &nb);

  /* Rescale diagonals */
  BLASC(scal)(&nb, SCS_BLAS_COMPLEX_CAST(&csqrt2_inv),
              SCS_BLAS_COMPLEX_CAST(w->cwork), &nb_p1);

  /* Repack into X */
  for (i = 0; i < n - 1; ++i) {
    X[i * (2 * n - i)] = w->cwork[i * (n + 1)][0];
    memcpy(&(X[i * (2 * n - i) + 1]), &(w->cwork[i * (n + 1) + 1]),
           2 * (n - i - 1) * sizeof(scs_float));
  }
  X[n * n - 1] = w->cwork[n * n - 1][0];
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

/* Box cone projection.
 *
 * proj_box_cone projects (t0, x) = tx onto
 *   K = { (t, s) : t l <= s <= t u, t >= 0 }
 * in the metric rho_t (t - t0)^2 + sum_j r_j (s_j - x_j)^2, where
 * rho_t = 1 / r_box[0] and r_j = 1 / r_box[j + 1] (all one without r_box).
 * For a fixed t the best s is x clipped to [t l, t u], so the projection is
 * the minimiser t* over t >= 0 of the convex piecewise quadratic
 *   f(t) = rho_t (t - t0)^2 / 2 + sum_j r_j dist(x_j, [t l_j, t u_j])^2 / 2.
 * Its derivative is
 *   g(t) = rho_t (t - t0) + sum_{j in U(t)} r_j (t u_j - x_j) u_j
 *                         + sum_{j in L(t)} r_j (t l_j - x_j) l_j,
 * where U(t) = { j : x_j > t u_j } and L(t) = { j : x_j < t l_j } are the rows
 * clipped at their upper and at their lower bound (never both, as l <= u).
 * g is continuous and strictly increasing (its slope is at least rho_t > 0),
 * so t* = 0 if g(0) >= 0 and otherwise t* is the root of g.
 *
 * Each of a row's three states (clipped above, clipped below, free) is cut out
 * by inequalities linear in t, so each holds on an interval of t. A row's
 * state therefore changes only at its kinks t = x_j / u_j and t = x_j / l_j,
 * and if a row has the same state at two points it has it at every point
 * between them. Between consecutive kinks, where U and L are fixed, g is the
 * line g(t) = H t - K with
 *   H = rho_t + sum_U r_j u_j^2 + sum_L r_j l_j^2,
 *   K = rho_t t0 + sum_U r_j x_j u_j + sum_L r_j x_j l_j,
 * so t* = K / H on the piece where g changes sign. Infinite bounds never clip
 * (t * inf = inf for t > 0; at t = 0, 0 * inf = NaN compares false, which
 * gives the closure of K there).
 *
 * The two methods below, A1 (box_t_sort) and A2 (box_t_newton), both end
 * with a bracket on which g changes sign, shown by g evaluated at its two
 * ends, and a single set of clipped rows that is right on the whole bracket,
 * so that g is one line there. They return the root of that line: A1 as
 * K / H, with H and K summed from scratch in one pass over the rows, A2 as
 * t - g(t) / H from a point t on that line close enough to the root that H's
 * rounding does not matter. The sums that set the result are compensated
 * (box_twosum), so the result is t* up to the rounding of the products in the
 * terms and of the final division, a bound that does not grow with the number
 * of rows:
 *   A1: |t - t*| <= ~3 eps (t* + (rho_t |t0| + sum |r_j x_j b_j|) / H),
 *   A2: |t - t*| <= ~3 eps (t* + sum |r_j (t b_j - x_j) b_j| / H),
 * the sums over the clipped rows (validated against an exact rational t*: at
 * most 4 eps max(t*, 1) max(1, cond), cond the A1 bound's sum over
 * H max(t*, 1)). There is no tolerance and no iteration limit that can stop
 * either short of the root. The two agree to that rounding, not bit for bit.
 * Every sign of g is taken from g in the residual form above, compensated,
 * with each row's state decided by the exact sign of t b - x (box_resid). A row
 * is thus never counted as clipped at a point where it is not, and a steep row
 * (large r b^2) contributes nothing at its own kink.
 *
 * Why exact: SCS 3.3's Newton loop stopped when |t - t_prev| or |g / H| fell
 * below about 1e-11 max(t, 1), and after 25 iterations. A step that small
 * does not show that t is near t*: when a steep row (large r b^2, e.g. a
 * bound of 1e20 left finite by normalize = 0) has its kink next to t, its
 * curvature makes the step tiny while t* lies beyond the kink. In
 * qafiro_tiny_qp without normalization (bounds of 1e20) it returned t of
 * about 1e-20 from a warm start near 0 where t* was 0.002 to 0.3: on 12 of
 * 2,931 calls with the 3.3 weights, from which ADMM recovered, and on 75,000
 * of 100,154 with the box t entry weighted like an equality row
 * (BOX_T_RY_FACTOR), where each such t was the next call's warm start and
 * ADMM cycled to its iteration limit. An inexact projection can stall ADMM
 * like that, so the methods below have no tolerance to meet: they stop only
 * on a sign change of g across an interval with no kink in it.
 *
 *   box_t_sort (A1): collect and sort the kinks, walk the pieces. O(m log m)
 *     per call (m = bsize - 1), no state carried between calls.
 *   box_t_newton (A2): Newton's method on g from the previous call's t, one
 *     O(m) pass per point, kept inside a bracket; stops when no kink lies
 *     between a point and the root of its line, or when two consecutive points
 *     on either side of the root have the same clipped sets. After
 *     BOX_NEWTON_MAX_PASSES passes, or a step that leaves the bracket, it hands
 *     the bracket to A1. Worst case O(m log m), typically one to three passes.
 *
 * proj_box_cone uses A2, with A1 as its fallback. Tradeoff: per call A2
 * costs about 1.7x the 3.3 loop on box cones dumped from real solves (the
 * compensated sums, the per-row states and the kink sweep; A1 alone costs
 * about 3.8x), which came to 2.5% more time per SCS iteration (geomean, -0.1%
 * to +5.5%) on the 10 benchmark problems with the largest box cones, where
 * the box projection is 2-20% of an iteration. In exchange no tolerance can
 * stop it short of t*. Validated against an exact rational t* on 5,849
 * instances (dumped SCS inputs and adversarial families: steep kinks next to
 * the root, bounds from 1e-6 to 1e12, ties, warm starts far off): within
 * 4 eps max(t*, 1) max(1, cond) on every one.
 */
/* A2: points evaluated (passes over the rows) before what is left of the
 * bracket goes to A1. Tradeoff: more passes let Newton finish more calls
 * without the O(m log m) sort, but a call that needs more than a few is one
 * where Newton is making little progress (slope of g dropping towards t*),
 * and A1 then finishes it in one sort. Inside SCS A2 took 1.9 passes per call
 * on average (53,514 calls, 14 problems); on the validation instances 8
 * passes ran out on 8 of 5,849. */
#define BOX_NEWTON_MAX_PASSES (8)
/* A2: a bound on the rounding of a Newton step's target tn, in eps times
 * (tn + sum |terms of g| / H): each term of g carries a few roundings, its
 * compensated sum one more, and the step and the division one each (see
 * box_t_newton). The next point is put this far past the predicted root.
 * Tradeoff: too small and that point can land on the same side of t* (no
 * sign change, another pass); too large and it can step past a nearby kink
 * (a state change, another pass). Either way it costs passes, not accuracy:
 * every stop is still a sign change of g across an interval with no kink. */
#define BOX_NEWTON_OVERSHOOT (16.)

typedef struct {
  scs_float t;  /* the kink */
  scs_int code; /* 2 j for the upper bound of row j, 2 j + 1 for its lower */
} BoxKink;

#ifdef SFLOAT
#define BOX_FMA fmaf
#define BOX_NEXT nextafterf
#define BOX_EPS FLT_EPSILON
#else
#define BOX_FMA fma
#define BOX_NEXT nextafter
#define BOX_EPS DBL_EPSILON
#endif
/* below this fraction of |x|, t b - x is formed with a single rounding (fma),
 * whose sign is exact. Above it the plain t * b - x cannot have the wrong
 * sign: its two roundings are at most a few eps |x|, far below |x| / 8.
 * Tradeoff: any fraction above a few eps would keep the sign exact; a larger
 * one sends more rows through the fma, which is slower on targets without a
 * hardware fma, a smaller one leaves less margin over the roundings. */
#define BOX_RESID_EXACT (0.125)

/* t b - x for a bound b: the residual of row j's term in g, < 0 exactly when
 * an upper bound b clips and > 0 exactly when a lower bound b clips. Near the
 * kink x / b the product and the difference cancel, and t * b - x can be wrong
 * by an ulp of x, sign included; so when |t b - x| < |x| / 8 it is formed with
 * one rounding (fma), whose sign is exact. Further out its two roundings cost
 * a few ulps of the result and cannot change its sign. An infinite bound gives
 * +-inf (t > 0) or NaN (t = 0), neither of which clips. */
static inline scs_float box_resid(scs_float t, scs_float b, scs_float x) {
  scs_float d = t * b - x;
  if (ABS(d) < BOX_RESID_EXACT * ABS(x))
    d = BOX_FMA(t, b, -x);
  return d;
}

/* (*s, *c) += v: Knuth's TwoSum, *s the rounded running sum and *c the sum of
 * its rounding errors, so that *s + *c is v's total to a few ulps whatever the
 * number of terms or the cancellation between them. Needs IEEE arithmetic (no
 * -ffast-math, which would drop the error terms). */
static inline void box_twosum(scs_float *s, scs_float *c, scs_float v) {
  scs_float t = *s + v, bp = t - *s;
  *c += (*s - (t - bp)) + (v - bp);
  *s = t;
}

/* H += r b^2, K += r x b, compensated (A1's last pass) */
static inline void box_hk(scs_float *h, scs_float *k, scs_float r, scs_float xj,
                          scs_float b) {
  box_twosum(&h[0], &h[1], r * b * b);
  box_twosum(&k[0], &k[1], r * xj * b);
}

/* One pass over the rows at t (A2): g(t) in the residual form, compensated
 * (box_twosum), the slope H of the piece at t, and *sg = sum |terms of g|, the
 * scale of g's rounding. With st, the state of every row is written there
 * (1 clipped above, -1 below, 0 free); with cmp also, the return value is the
 * number of rows whose state differs from the one st held (else 0). */
static scs_int box_pass(const scs_float *x, const scs_float *bl,
                        const scs_float *bu, scs_int n, scs_float t0,
                        scs_float rho_t, const scs_float *rinv, scs_float t,
                        signed char *st, scs_int cmp, scs_float *g,
                        scs_float *h, scs_float *sg) {
  scs_float gt[2] = {rho_t * (t - t0), 0.}, ht = rho_t, sgt = ABS(gt[0]), r,
            b = 0., d, v;
  scs_int j, nchg = 0;
  signed char sj;
  for (j = 0; j < n; ++j) {
    sj = 0;
    d = box_resid(t, bu[j], x[j]);
    if (d < 0.) {
      b = bu[j];
      sj = 1;
    } else {
      d = box_resid(t, bl[j], x[j]);
      if (d > 0.) {
        b = bl[j];
        sj = -1;
      }
    }
    if (sj) {
      r = rinv ? rinv[j] : 1.0;
      v = r * d * b;
      box_twosum(&gt[0], &gt[1], v);
      sgt += ABS(v);
      ht += r * b * b;
    }
    if (st) {
      if (cmp)
        nchg += (st[j] != sj);
      st[j] = sj;
    }
  }
  *g = gt[0] + gt[1];
  *h = ht;
  if (sg)
    *sg = sgt;
  return nchg;
}

/* 1 if some row has a kink within s of t, at which its state could change.
 * The kink of bound b is at distance |t b - x| / |b| from t; the test is
 * |t b - x| <= 2 s |b|, the factor 2 covering the rounding of t b - x (at most
 * eps |t b| when not fused, and s >= 16 eps t where it is used). An infinite
 * bound (t b - x infinite, or NaN at t = 0) has no kink; a zero one has none
 * either unless x = 0, when it is reported, which is only conservative. No
 * branch and no state carried between rows except the flag, so the loop
 * vectorises: it costs a fraction of a box_pass. */
static scs_int box_kink_near(const scs_float *x, const scs_float *bl,
                             const scs_float *bu, scs_int n, scs_float t,
                             scs_float s) {
  scs_float s2 = 2. * s, du, dl;
  scs_int j, found = 0;
  for (j = 0; j < n; ++j) {
    du = ABS(t * bu[j] - x[j]);
    dl = ABS(t * bl[j] - x[j]);
    found |= (du <= s2 * ABS(bu[j])) & (du < INFINITY);
    found |= (dl <= s2 * ABS(bl[j])) & (dl < INFINITY);
  }
  return found;
}

/* g(t) in the residual form, with each row's state at t (for the checks) */
static inline scs_float box_g_term(scs_float t, scs_float xj, scs_float l,
                                   scs_float u, scs_float r) {
  scs_float d = box_resid(t, u, xj);
  if (d < 0.)
    return r * d * u;
  d = box_resid(t, l, xj);
  if (d > 0.)
    return r * d * l;
  return 0.;
}

/* The last pass of A1, on the piece (a, b) that holds m: H and K over the
 * clipped sets at m, and g at a and at b, each with the rows' states there
 * (b = inf: gb is not formed). */
static void box_piece(const scs_float *x, const scs_float *bl,
                      const scs_float *bu, scs_int n, scs_float t0,
                      scs_float rho_t, const scs_float *rinv, scs_float a,
                      scs_float m, scs_float b, scs_float *h, scs_float *k,
                      scs_float *ga, scs_float *gb) {
  scs_float ht[2] = {rho_t, 0.}, kt[2] = {rho_t * t0, 0.},
            gat[2] = {rho_t * (a - t0), 0.}, gbt[2] = {rho_t * (b - t0), 0.}, r;
  scs_int j, fb = b < INFINITY;
  for (j = 0; j < n; ++j) {
    r = rinv ? rinv[j] : 1.0;
    if (box_resid(m, bu[j], x[j]) < 0.)
      box_hk(ht, kt, r, x[j], bu[j]);
    else if (box_resid(m, bl[j], x[j]) > 0.)
      box_hk(ht, kt, r, x[j], bl[j]);
    box_twosum(&gat[0], &gat[1], box_g_term(a, x[j], bl[j], bu[j], r));
    if (fb)
      box_twosum(&gbt[0], &gbt[1], box_g_term(b, x[j], bl[j], bu[j], r));
  }
  *h = ht[0] + ht[1];
  *k = kt[0] + kt[1];
  *ga = gat[0] + gat[1];
  *gb = fb ? gbt[0] + gbt[1] : INFINITY;
}

/* Neumaier's compensated sum: (*s + *c) += v. The walk over the kinks adds and
 * removes rows from H and K; this keeps a large term that leaves from
 * swamping the small ones that stay. */
static inline void box_csum(scs_float *s, scs_float *c, scs_float v) {
  scs_float t = *s + v;
  if (ABS(*s) >= ABS(v))
    *c += (*s - t) + v;
  else
    *c += (v - t) + *s;
  *s = t;
}

/* Introsort of the kinks by t: quicksort (median of three, Hoare partition),
 * insertion sort below 16 entries, heapsort past 2 log2(n) levels, so
 * O(n log n) in the worst case. Deterministic. */
static void box_kink_isort(BoxKink *a, scs_int n) {
  scs_int i, j;
  BoxKink v;
  for (i = 1; i < n; ++i) {
    v = a[i];
    for (j = i; j > 0 && a[j - 1].t > v.t; --j)
      a[j] = a[j - 1];
    a[j] = v;
  }
}

static void box_kink_sift(BoxKink *a, scs_int i, scs_int n) {
  BoxKink v = a[i];
  scs_int c;
  while ((c = 2 * i + 1) < n) {
    if (c + 1 < n && a[c + 1].t > a[c].t)
      ++c;
    if (!(a[c].t > v.t))
      break;
    a[i] = a[c];
    i = c;
  }
  a[i] = v;
}

static void box_kink_qsort(BoxKink *a, scs_int n, scs_int depth) {
  BoxKink v;
  scs_float p;
  scs_int i, j, m;
  while (n > 16) {
    if (depth-- <= 0) { /* heapsort */
      for (i = n / 2 - 1; i >= 0; --i)
        box_kink_sift(a, i, n);
      for (i = n - 1; i > 0; --i) {
        v = a[0];
        a[0] = a[i];
        a[i] = v;
        box_kink_sift(a, 0, i);
      }
      return;
    }
    m = n / 2;
    if (a[m].t < a[0].t) {
      v = a[m];
      a[m] = a[0];
      a[0] = v;
    }
    if (a[n - 1].t < a[0].t) {
      v = a[n - 1];
      a[n - 1] = a[0];
      a[0] = v;
    }
    if (a[n - 1].t < a[m].t) {
      v = a[n - 1];
      a[n - 1] = a[m];
      a[m] = v;
    }
    p = a[m].t;
    i = -1;
    j = n;
    for (;;) {
      do
        ++i;
      while (a[i].t < p);
      do
        --j;
      while (p < a[j].t);
      if (i >= j)
        break;
      v = a[i];
      a[i] = a[j];
      a[j] = v;
    }
    /* a[0 .. j] <= p <= a[j + 1 .. n - 1], both parts non-empty: recurse on
     * the shorter part, loop on the longer */
    if (j + 1 < n - j - 1) {
      box_kink_qsort(a, j + 1, depth);
      a += j + 1;
      n -= j + 1;
    } else {
      box_kink_qsort(a + j + 1, n - j - 1, depth);
      n = j + 1;
    }
  }
  box_kink_isort(a, n);
}

static void box_kink_sort(BoxKink *a, scs_int n) {
  scs_int depth = 0, m = n;
  while (m > 1) {
    m >>= 1;
    depth += 2;
  }
  box_kink_qsort(a, n, depth);
}

/* Bound b of row j clips on the interval (ca, cb) of t > 0 (ca in {0, kink},
 * cb in {kink, inf}); returns 0 if it never clips. upper: it clips where
 * x > t b, else where x < t b. */
static scs_int box_clip_interval(scs_float xj, scs_float b, scs_int upper,
                                 scs_float *ca, scs_float *cb) {
  if (!(ABS(b) < INFINITY))
    return 0;
  if (upper ? b > 0. : b < 0.) { /* clips for t below the kink */
    if (upper ? !(xj > 0.) : !(xj < 0.))
      return 0;
    *ca = 0.;
    *cb = xj / b;
  } else if (b != 0.) { /* clips for t above the kink */
    *ca = (upper ? xj >= 0. : xj <= 0.) ? 0. : xj / b;
    *cb = INFINITY;
  } else { /* b = 0: clips for every t > 0 or never */
    if (upper ? !(xj > 0.) : !(xj < 0.))
      return 0;
    *ca = 0.;
    *cb = INFINITY;
  }
  return 1;
}

/* First index k in [k0, k1) with g(kk[k].t) >= 0 (k1 if none), by bisection
 * with g in the residual form */
static scs_int box_kink_search(const scs_float *x, const scs_float *bl,
                               const scs_float *bu, scs_int n, scs_float t0,
                               scs_float rho_t, const scs_float *rinv,
                               const BoxKink *kk, scs_int k0, scs_int k1) {
  scs_float g, h;
  scs_int mid;
  while (k0 < k1) {
    mid = k0 + (k1 - k0) / 2;
    box_pass(x, bl, bu, n, t0, rho_t, rinv, kk[mid].t, SCS_NULL, 0, &g, &h,
             SCS_NULL);
    if (g >= 0.)
      k1 = mid;
    else
      k0 = mid + 1;
  }
  return k0;
}

/* A1: t* in [lo, hi], given hi = inf or g(hi) >= 0, and g(lo) < 0 if lo_neg
 * (else lo = 0 and the sign of g(0) is not known yet). kk: workspace of
 * 2 n kinks.
 * 1. One pass: the kinks in (lo, hi), and H and K on the first piece (lo, .).
 * 2. lo = 0 and not lo_neg: g(0) = -K, so t* = 0 if K <= 0.
 * 3. Sort the kinks and walk them in order, updating H and K (compensated) as
 *    rows enter and leave, until g = H c - K >= 0 at a kink c. At c, the rows
 *    whose kink c is are left out of g: their term is zero there. t* is then
 *    on the piece (a, b) that ends at c.
 * 4. One pass (box_piece): H and K from scratch over the clipped sets at the
 *    midpoint m of (a, b), and g at a and at b with the rows' own states
 *    there. If those say t* is outside [a, b] (the walk's running sums can
 *    misjudge a sign only within their rounding), the piece is found again by
 *    bisection over the sorted kinks with g in the residual form, O(log m)
 *    passes, and step 4 is repeated once on it without the check.
 * 5. Return K / H clamped to [a, b]. The clamp matters only when a kink lies
 *    within an ulp inside a or b, so that the line of m does not hold right
 *    up to that end: g changes sign on [a, b], so t* is then within that ulp of
 *    the end it is clamped to. If no double lies strictly inside (a, b), t* is
 *    within an ulp of b, which is returned.
 * Termination: each step is a pass over the rows, a sort, a walk over a finite
 * list or a bisection over it, and step 4 runs at most twice. */
static scs_float box_t_sort(const scs_float *x, const scs_float *bl,
                            const scs_float *bu, scs_int n, scs_float t0,
                            scs_float rho_t, const scs_float *rinv,
                            scs_float lo, scs_int lo_neg, scs_float hi,
                            BoxKink *kk) {
  scs_float hs = rho_t, hc = 0., ks = rho_t * t0, kc = 0., r, ca, cb, c, a,
            b, m, h, k, ga, gb, bd, t;
  scs_int j, ne = 0, i, i2, code, lower, repaired;
  for (j = 0; j < n; ++j) {
    for (lower = 0; lower < 2; ++lower) {
      bd = lower ? bl[j] : bu[j];
      if (!box_clip_interval(x[j], bd, !lower, &ca, &cb))
        continue;
      if (ca <= lo && lo < cb) {
        r = rinv ? rinv[j] : 1.0;
        box_csum(&hs, &hc, r * bd * bd);
        box_csum(&ks, &kc, r * x[j] * bd);
      }
      code = 2 * j + lower;
      if (ca > lo && ca < hi) {
        kk[ne].t = ca;
        kk[ne++].code = code;
      }
      if (cb > lo && cb < hi) {
        kk[ne].t = cb;
        kk[ne++].code = code;
      }
    }
  }
  if (!lo_neg && !(ks + kc > 0.)) { /* g(0) = -K >= 0 */
    return 0.;
  }
  box_kink_sort(kk, ne);
  i = 0;
  while (i < ne) {
    c = kk[i].t;
    for (i2 = i; i2 < ne && kk[i2].t == c; ++i2) {
      code = kk[i2].code;
      j = code >> 1;
      lower = code & 1;
      bd = lower ? bl[j] : bu[j];
      if (!(lower ? bd > 0. : bd < 0.)) { /* leaves at c */
        r = rinv ? rinv[j] : 1.0;
        box_csum(&hs, &hc, -(r * bd * bd));
        box_csum(&ks, &kc, -(r * x[j] * bd));
      }
    }
    if ((hs + hc) * c - (ks + kc) >= 0.)
      break;
    for (; i < i2; ++i) {
      code = kk[i].code;
      j = code >> 1;
      lower = code & 1;
      bd = lower ? bl[j] : bu[j];
      if (lower ? bd > 0. : bd < 0.) { /* enters at c */
        r = rinv ? rinv[j] : 1.0;
        box_csum(&hs, &hc, r * bd * bd);
        box_csum(&ks, &kc, r * x[j] * bd);
      }
    }
  }
  /* t* on the piece (a, b) = (kk[i - 1].t, kk[i].t), or ending at lo / hi */
  for (repaired = 0;; repaired = 1) {
    a = i > 0 ? kk[i - 1].t : lo;
    b = i < ne ? kk[i].t : hi;
    m = b < INFINITY ? a + (b - a) / 2 : (a > 0. ? 2 * a : 1.);
    if (!(m > a && m < b)) /* no double strictly inside */
      return b;
    box_piece(x, bl, bu, n, t0, rho_t, rinv, a, m, b, &h, &k, &ga, &gb);
    if (!repaired && b < hi && gb < 0.) { /* t* > b */
      i = box_kink_search(x, bl, bu, n, t0, rho_t, rinv, kk, i + 1, ne);
    } else if (!repaired && i > 0 && ga >= 0.) { /* t* <= a */
      i = box_kink_search(x, bl, bu, n, t0, rho_t, rinv, kk, 0, i - 1);
    } else {
      t = k / h;
      t = MIN(MAX(t, a), b);
      return MAX(t, 0.);
    }
  }
}

/* A2: Newton's method on g from the warm start t_wm (the previous call's t*),
 * kept inside a bracket [lo, hi] on t*. Each point t costs one pass over the
 * rows (box_pass): g(t) in the residual form (compensated), the slope H of
 * t's piece, and the rows' states, compared with those at the previous point.
 * The root of the line of t's piece is tn = t - g / H (0 if that is
 * negative), up to a rounding bound err = BOX_NEWTON_OVERSHOOT eps
 * (tn + sum |terms of g| / H).
 * - Stop, one point: no kink within |tn - t| + err of t (box_kink_near, a
 *   branch-free sweep over the rows), so g is that line up to tn and tn is t*
 *   (if clamped at 0, g(0) >= 0 and t* = 0). Also required:
 *   n |tn - t| <= max(tn, 1), so that the rounding of H, a plain sum (at most
 *   n eps relative), moves tn by under an ulp of max(t*, 1). This is the
 *   usual end when the previous call's t* is on the right piece and no kink
 *   is nearer to t than the step: one pass and one sweep.
 * - Stop, two points: g < 0 at one of two consecutive points and g >= 0 at the
 *   other, with the same state in every row. Every row then keeps its state
 *   between them (see above), so g is one line there and t* = tn from the
 *   second point, which lies within the overshoot below of t*: H's rounding
 *   does not matter there.
 * - Adjacent doubles lo < hi with g(lo) < 0 <= g(hi): t* is within an ulp of
 *   hi, which is returned (as A1 does).
 * - g(0) >= 0 at t = 0: t* = 0.
 * - Otherwise the next point is tn moved on by err past the root it predicts
 *   (towards hi when g < 0, towards lo when g >= 0). The step lands within its
 *   rounding of that root, where g has either sign; the overshoot puts the
 *   point on the far side, so that when the root is on the current piece the
 *   next pass sees g change sign with the same states and stops. A step below
 *   the rounding of t still moves the point, so the root is confirmed by a
 *   sign change: a step that small does not show that the root is near t when
 *   a steep row's kink lies within that rounding of t (the trap of the
 *   original loop's stopping test). A point that would fall outside the
 *   bracket is taken without the overshoot.
 * - Every point narrows [lo, hi] by the sign of g there. A step that does not
 *   land strictly inside the bracket (Newton overshoots where the slope of g
 *   drops towards t*), or BOX_NEWTON_MAX_PASSES passes, hands the bracket to
 *   A1, which finds t* in it.
 * Termination: at most BOX_NEWTON_MAX_PASSES passes, then A1, which
 * terminates (see box_t_sort). The states of the previous point are kept in
 * A1's workspace, which A1 only uses after the last pass. */
static scs_float box_t_newton(const scs_float *x, const scs_float *bl,
                              const scs_float *bu, scs_int n, scs_float t0,
                              scs_float rho_t, const scs_float *rinv,
                              scs_float t_wm, BoxKink *kk) {
  scs_float t = t_wm, g, h, sg, tn, to, err, gp = 0., lo = 0., hi = INFINITY;
  scs_int pass, lo_neg = 0, nchg;
  signed char *st = (signed char *)kk;
  if (!(t >= 0. && t < INFINITY))
    t = 1.;
  for (pass = 1;; ++pass) {
    nchg = box_pass(x, bl, bu, n, t0, rho_t, rinv, t, st, pass > 1, &g, &h,
                    &sg);
    /* the root of the line of t's piece (0 if that is negative), and a bound
     * on its rounding */
    tn = MAX(t - g / h, 0.);
    err = BOX_NEWTON_OVERSHOOT * BOX_EPS * (tn + sg / h);
    if ((pass > 1 && nchg == 0 && (g >= 0.) != (gp >= 0.)) ||
        ((scs_float)n * ABS(tn - t) <= MAX(tn, 1.) &&
         !box_kink_near(x, bl, bu, n, t, ABS(tn - t) + err))) {
      /* g changed sign since the previous point with no state change, which
       * puts the root within the overshoot of t; or no kink lies within the
       * step and its rounding, so g is t's line up to the root (g(0) >= 0 if
       * clamped), and the step is short enough that H's rounding (H a plain
       * sum, at most n eps relative) moves the root by under an ulp of
       * max(t*, 1) */
      return tn;
    }
    if (g >= 0.) {
      if (t <= 0.) {
        return 0.;
      }
      hi = t;
    } else {
      lo = t;
      lo_neg = 1;
    }
    if (lo_neg && hi < INFINITY && BOX_NEXT(lo, hi) == hi) {
      return hi;
    }
    /* the Newton step, moved on past the root it predicts (towards hi when
     * g < 0, towards lo when g >= 0) */
    to = MAX(g < 0. ? tn + err : tn - err, 0.);
    if (to < hi && (to > lo || (to == 0. && !lo_neg)))
      tn = to;
    if (!(tn < hi && (tn > lo || (tn == 0. && !lo_neg)))) {
      break;
    }
    if (pass >= BOX_NEWTON_MAX_PASSES) {
      break;
    }
    t = tn;
    gp = g;
  }
  return box_t_sort(x, bl, bu, n, t0, rho_t, rinv, lo, lo_neg, hi, kk);
}

/* Project onto { (t, s) | t * l <= s <= t * u, t >= 0 } (see above).
   tx = [t; s], total length = bsize, under Euclidean metric 1/r_box.
*/
static scs_float proj_box_cone(scs_float *tx, const scs_float *bl,
                               const scs_float *bu, scs_int bsize,
                               scs_float t_wm, scs_float *r_box,
                               scs_float *r_box_inv, BoxKink *kk) {
  scs_float *x = &(tx[1]);
  scs_float t, rho_t = 1.0;
  const scs_float *rinv = SCS_NULL;
  scs_int j;

  if (bsize == 1) {
    tx[0] = MAX(tx[0], 0.0);
    return tx[0];
  }
  if (r_box) {
    rho_t = 1.0 / r_box[0];
    /* Precompute reciprocals once, used across all passes */
    if (r_box_inv) {
      for (j = 0; j < bsize - 1; j++) {
        r_box_inv[j] = 1.0 / r_box[j + 1];
      }
      rinv = r_box_inv;
    }
  }

  t = box_t_newton(x, bl, bu, bsize - 1, tx[0], rho_t, rinv, t_wm, kk);

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
  scs_float v1, s, alpha;
  if (q <= 0)
    return;
  if (q == 1) {
    x[0] = MAX(x[0], 0.);
    return;
  }

  v1 = x[0];
  /* Fast paths for the two most common small SOC sizes avoid BLAS call
   * overhead (function pointer dispatch + Fortran ABI arguments). */
  if (q == 2) {
    s = ABS(x[1]);
  } else if (q == 3) {
    s = SQRTF(x[1] * x[1] + x[2] * x[2]);
  } else {
    s = SCS(norm_2)(&(x[1]), q - 1);
  }
  alpha = (s + v1) / 2.0;

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
  /* Positive root of x^2 - xh * x - a * (rh - r) * r = 0, in a form that
   * avoids cancellation when xh < 0. */
  scs_float d = 4 * a * (rh - r) * r, s = SQRTF(xh * xh + d);
  return xh >= 0 ? 0.5 * (xh + s) : 0.5 * d / (s - xh);
}

static void proj_power_cone(scs_float *v, scs_float a) {
  scs_float xh = v[0], yh = v[1], rh = ABS(v[2]);
  scs_float x = 0.0, y = 0.0, r, lo, hi;
  scs_float f_lo, f_hi; /* f at the bracket ends, for regula falsi */
  scs_float phi0; /* phi(0) = phi(rh) = max(xh, 0)^a * max(yh, 0)^(1-a) */
  /* Projection onto a cone is positively homogeneous, so work on v scaled
   * to unit inf-norm; this makes the tolerances below relative to |v|. */
  scs_float nrm = MAX(MAX(ABS(xh), ABS(yh)), rh);
  scs_int i, side = 0;

  if (nrm <= 0)
    return;
  xh /= nrm;
  yh /= nrm;
  rh /= nrm;

  if (xh >= 0 && yh >= 0)
    phi0 = POWF(xh, a) * POWF(yh, (1 - a));
  else if (a > 0 && a < 1)
    phi0 = 0;
  else
    phi0 = POWF(MAX(xh, 0), a) * POWF(MAX(yh, 0), (1 - a));

  /* Check v membership in K_a */
  if (xh >= 0 && yh >= 0 && POW_CONE_TOL + phi0 >= rh)
    return;

  /* Check -v membership in Polar K_a^* */
  if (xh <= 0 && yh <= 0 &&
      POW_CONE_TOL + POWF(-xh, a) * POWF(-yh, 1 - a) >=
          rh * POWF(a, a) * POWF(1 - a, 1 - a)) {
    v[0] = v[1] = v[2] = 0;
    return;
  }

  /* Find the root r in (0, rh) of f(r) = x(r)^a * y(r)^(1-a) - r. The first
   * term is concave in r, so f > 0 left of the root and f < 0 right of it,
   * which gives a bracket [lo, hi] that safeguards the Newton steps. f_lo
   * and f_hi are f at the bracket ends, f(0) = phi(0) and f(rh) = phi(rh) - rh
   * to start, for the regula falsi fallback. */
  lo = 0;
  hi = rh;
  f_lo = phi0;
  f_hi = phi0 - rh;
  r = rh / 2;
  /* With exactly one negative component, for small r x(r) ~ a * rh * r / |xh|
   * (or the same for y), and phi(r) = r then gives the root directly. This is
   * the common case for ADMM iterates, and its root is often tiny, which
   * Newton from rh / 2 reaches only slowly. */
  if (xh < 0 && yh > 0) {
    scs_float r0 = POWF(a * rh / (-xh), a / (1 - a)) * yh;
    if (r0 > 0 && r0 < rh / 2)
      r = r0;
  } else if (yh < 0 && xh > 0) {
    scs_float r0 = POWF((1 - a) * rh / (-yh), (1 - a) / a) * xh;
    if (r0 > 0 && r0 < rh / 2)
      r = r0;
  }
  /* Points just outside the cone (e.g. converging iterates) have their root
   * near rh; starting there avoids Newton overshooting from the left into
   * bisection. Only do this when the point is close to the cone,
   * phi(rh) >= rh / 2, and the slopes dx/dr = -a * rh / xh and
   * dy/dr = -(1 - a) * rh / yh at rh are moderate: otherwise the first step
   * is too small to represent and looks converged. */
  if (xh > 0 && yh > 0 && phi0 >= rh / 2 &&
      a * rh <= POW_CONE_START_SLOPE * xh &&
      (1 - a) * rh <= POW_CONE_START_SLOPE * yh)
    r = rh;
  for (i = 0; i < POW_CONE_MAX_ITERS; ++i) {
    scs_float f, fp, dxdr, dydr, xa, y1a, r_new, r_g, sens;
    x = pow_calc_x(r, xh, rh, a);
    y = pow_calc_x(r, yh, rh, 1 - a);

    /* Cache POWF(x,a) and POWF(y,1-a): both are needed for f and fp,
     * so computing them once saves two POWF calls per Newton iteration. */
    xa  = POWF(x, a);
    y1a = POWF(y, (1 - a));
    f = xa * y1a - r;
    /* Illinois: halve the f value kept at an end that was not just moved
     * twice in a row, so regula falsi does not stagnate. */
    if (f > 0) {
      lo = r;
      f_lo = f;
      if (side > 0)
        f_hi *= 0.5;
      side = 1;
    } else {
      hi = r;
      f_hi = f;
      if (side < 0)
        f_lo *= 0.5;
      side = -1;
    }

    dxdr = a * (rh - 2 * r) / (2 * x - xh);
    dydr = (1 - a) * (rh - 2 * r) / (2 * y - yh);
    fp = xa * y1a * (a * dxdr / x + (1 - a) * dydr / y) - 1;

    r_new = r - f / fp;
    if (f < 0) {
      /* Near the polar boundary f is nearly flat around a double root and
       * Newton on f converges only linearly; Newton on f(r) / r, which
       * removes the trivial root at 0, converges much faster there. It is
       * also taken elsewhere: restricting it saved a few projection
       * evaluations on some inputs but made whole solves worse. */
      r_g = r - r * f / (r * fp - f);
      if (r_g > lo && r_g < hi)
        r_new = r_g;
    }
    /* Stop when the step is small relative to r and moves x, y and r by at
     * most POW_CONE_TOL (relative to |v| after normalization). Stop on the
     * step rather than on |f|: a small |f| does not imply r is accurate when
     * f is flat, and near the apex x can be far more sensitive to r than r
     * is. r is an end of the bracket here, so a converged step may land
     * exactly on it. The relative test cannot be met when the root is at or
     * near 0, so also stop once the bracket is that narrow, or when it cannot
     * shrink any further in the working precision. */
    sens = MAX(1, MAX(ABS(dxdr), ABS(dydr)));
    if (ABS(r_new - r) <= POW_CONE_TOL * r &&
        ABS(r_new - r) * sens <= POW_CONE_TOL && r_new >= lo && r_new <= hi)
      break;
    /* If the Newton step leaves the bracket (from the left of a concave f it
     * overshoots), fall back to regula falsi, then to bisection, which also
     * catches NaN. Bisection alone is very slow when the root is within a few
     * ulps of rh, e.g. a point just outside the cone with a tiny component. */
    if (!(r_new > lo && r_new < hi))
      r_new = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
    if (!(r_new > lo && r_new < hi))
      r_new = 0.5 * (lo + hi);
    if ((hi - lo) * sens <= POW_CONE_TOL || hi - lo <= 4 * POW_CONE_EPS * hi)
      break;
    r = r_new;
  }
  v[0] = x * nrm;
  v[1] = y * nrm;
  v[2] = ((v[2] < 0) ? -r : r) * nrm;
}

/* ===================== PSD Blocks, all Threads ======================== */

/* Project the blocks of lists t, t + dt, t + 2 dt, ..., each list with its own
 * workspace: all of them in order for t = 0 and dt = 1, the share of thread t
 * in a team of dt otherwise. x and rho as in proj_psd_cones. */
static scs_int proj_psd_lists(scs_float *x, ScsConeWork *c,
                              const scs_float *rho, scs_int t, scs_int dt) {
  const ScsCone *k = c->k;
  scs_int j, b, st, status = 0;
  for (; t < c->psd_nthreads; t += dt) {
    ScsPsdWork *w = SCS_NULL; /* no workspace without LAPACK */
#ifdef USE_LAPACK
    w = &(c->psd_tw[t]);
#endif
    for (j = c->psd_tptr[t]; j < c->psd_tptr[t + 1]; ++j) {
      const ScsPsdLapackSizes *sz = SCS_NULL; /* none without LAPACK */
      b = c->psd_blk[j];
#ifdef USE_LAPACK
      sz = &(c->psd_lsz[b]);
#endif
      if (b < k->ssize) {
        st = proj_semi_definite_cone(&(x[c->psd_off[b]]), k->s[b], w, sz,
                                     rho ? &(rho[c->psd_off[b]]) : SCS_NULL);
      } else {
        st = proj_complex_semi_definite_cone(&(x[c->psd_off[b]]),
                                             k->cs[b - k->ssize], w, sz);
      }
      if (st < 0)
        status = st;
    }
  }
  return status;
}

/* Project the PSD blocks x (from the first PSD entry); rho is the matching
 * part of r_y or NULL. */
static scs_int proj_psd_cones(scs_float *x, ScsConeWork *c,
                              const scs_float *rho) {
#ifdef PSD_PARALLEL
  /* the threads of the split, or fewer if the OpenMP thread count has been
   * lowered since setup */
  int nt = omp_get_max_threads();
  if ((scs_int)nt > c->psd_nthreads)
    nt = (int)c->psd_nthreads;
  if (nt > 1) {
    scs_int status = 0;
#pragma omp parallel num_threads(nt)
    {
      /* thread t of a team of dt, which may be smaller than asked for (e.g.
       * nested), projects lists t, t + dt, ...: every list is covered. The
       * BLAS has to run on one thread in here (see the comment on the option
       * above). SCS can only see to that for MKL: one MKL thread meanwhile,
       * then back to the setting before. */
      int mkl_nt = mkl_set_local_threads(1);
      scs_int st = proj_psd_lists(x, c, rho, omp_get_thread_num(),
                                  omp_get_num_threads());
      mkl_set_local_threads(mkl_nt);
      if (st < 0) {
#pragma omp critical(scs_psd_status)
        status = st;
      }
    }
    return status;
  }
#endif
  /* no parallel region: the BLAS is free to use its own threads, as in a
   * build without USE_PARALLEL_PSD */
  return proj_psd_lists(x, c, rho, 0, 1);
}

/* ================ Main Cone Projection Dispatch ====================== */

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
                                        c->box_t_warm_start, r_box,
                                        c->r_box_inv, (BoxKink *)c->box_kinks);
    count += k->bsize;
  }

  /* 4. SOC */
  if (k->qsize) {
    for (i = 0; i < k->qsize; ++i) {
      proj_soc(&(x[count]), k->q[i]);
      count += k->q[i];
    }
  }

  /* 5-6. PSD and Complex PSD Cones */
  if (c->psd_nthreads > 0) {
#ifdef USE_SPECTRAL_CONES
    SPECTRAL_TIMING(SCS(tic)(&spec_mat_proj_timer);)
#endif
    status = proj_psd_cones(&(x[count]), c, r_y ? &(r_y[count]) : SCS_NULL);
#ifdef USE_SPECTRAL_CONES
    SPECTRAL_TIMING(c->tot_time_mat_cone_proj +=
                    SCS(tocq)(&spec_mat_proj_timer);)
#endif
    if (status < 0)
      return status;
    count += c->psd_len;
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

/* =========================== Public API ============================== */

ScsConeWork *SCS(init_cone)(ScsCone *k, scs_int m) {
  ScsConeWork *c = (ScsConeWork *)scs_calloc(1, sizeof(ScsConeWork));
  if (!c)
    return SCS_NULL;

  c->k = k;
  c->m = m;
  c->scaled_cones = 0;

  set_cone_boundaries(k, c);
  c->s = (scs_float *)scs_calloc(m, sizeof(scs_float));

  /* Preallocate workspace for box cone reciprocal metric */
  if (k->bsize > 1) {
    c->r_box_inv =
        (scs_float *)scs_calloc(k->bsize - 1, sizeof(scs_float));
    /* kinks of g for the box projection's sort (A1), at most two per row;
     * A2 keeps the rows' states in the same space */
    c->box_kinks = scs_calloc(2 * (k->bsize - 1), sizeof(BoxKink));
    if (!c->r_box_inv || !c->box_kinks) {
      SCS(finish_cone)(c);
      return SCS_NULL;
    }
  }

  /* PSD blocks: offsets and threads (the workspaces are set up below) */
  if ((k->ssize && k->s) || (k->cssize && k->cs)) {
    if (set_up_psd_schedule(c, k) < 0) {
      SCS(finish_cone)(c);
      return SCS_NULL;
    }
  }

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

/* Discard inner-solver state that was cached under a previous metric.
 *
 * The standard cones project in closed form, so their projections are a
 * function of the point and the current metric alone. The spectral cones
 * are not: proj_logdet_cone runs a Newton/IPM inner solve and warm-starts
 * Newton from the projection it produced on the previous outer iteration.
 * That saved iterate belongs to the metric that produced it -- once diag_r
 * changes, the point being projected jumps, and the stale iterate is no
 * longer a good (nor necessarily a domain-interior) starting point. So
 * callers must invalidate it whenever the metric changes, and whenever a
 * new solve begins.
 *
 * Clearing the flags is a complete invalidation: saved_log_projs is read
 * only where the matching flag is set, log_cone_Newton writes every
 * component of its output on the cold path, and the IPM fallback always
 * starts from the all-ones point, so it carries nothing across iterations.
 *
 * Dropping a warm start only costs inner iterations, never correctness,
 * which is why this is safe to call unconditionally.
 */
void SCS(reset_cone_cache)(ScsConeWork *c) {
#ifdef USE_SPECTRAL_CONES
  if (c && c->k && c->log_cone_warmstarts && c->k->dsize > 0) {
    memset(c->log_cone_warmstarts, 0, c->k->dsize * sizeof(bool));
  }
#else
  (void)c;
#endif
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

  /* x -> - Rx; hoist the r_y != NULL check outside the loop */
  if (r_y) {
    for (i = 0; i < c->m; ++i) {
      x[i] *= -r_y[i];
    }
  } else {
    for (i = 0; i < c->m; ++i) {
      x[i] = -x[i];
    }
  }

  /* Project -x onto cone, x -> \Pi_{C^*}^{R^{-1}}(-x) under r_y metric */
  status = proj_cone(x, k, c, scal ? 1 : 0, r_y);

  /* Return x + R^{-1} \Pi_{C^*}^{R^{-1}} ( -x ); hoist r_y check */
  if (r_y) {
    for (i = 0; i < c->m; ++i) {
      x[i] = x[i] / r_y[i] + c->s[i];
    }
  } else {
    for (i = 0; i < c->m; ++i) {
      x[i] += c->s[i];
    }
  }

  return status;
}
