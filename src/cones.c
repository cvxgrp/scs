#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "scs_blas.h" /* contains BLAS(X) macros and type info */
#include "util.h"

#define CONE_TOL (1e-9)
#define CONE_THRESH (1e-8)
#define EXP_CONE_MAX_ITERS (100)
#define BOX_CONE_MAX_ITERS (25)
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
blas_int BLAS(syrk)(const char *uplo, const char *trans, const blas_int *n,
                    const blas_int *k, const scs_float *alpha,
                    const scs_float *a, const blas_int *lda,
                    const scs_float *beta, scs_float *c, const blas_int *ldc);
void BLAS(scal)(const blas_int *n, const scs_float *sa, scs_float *sx,
                const blas_int *incx);

#ifdef __cplusplus
}
#endif

#endif

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
    if (k->p)
      scs_free(k->p);
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
  /* copy power cone */
  if (src->psize > 0) {
    dest->p = (scs_float *)scs_calloc(src->psize, sizeof(scs_float));
    memcpy(dest->p, src->p, src->psize * sizeof(scs_float));
  } else {
    dest->p = SCS_NULL;
  }
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

/*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size strictly
 * larger than 1, boundaries malloc-ed here so should be freed.
 */
void set_cone_boundaries(const ScsCone *k, ScsConeWork *c) {
  scs_int i, s_cone_sz, count = 0;
  scs_int cone_boundaries_len =
      1 + k->qsize + k->ssize + k->ed + k->ep + k->psize;
  scs_int *b = (scs_int *)scs_calloc(cone_boundaries_len, sizeof(scs_int));
  /* cones that can be scaled independently */
  b[count] = k->z + k->l + k->bsize;
  count += 1; /* started at 0 now move to first entry */
  for (i = 0; i < k->qsize; ++i) {
    b[count + i] = k->q[i];
  }
  count += k->qsize;
  for (i = 0; i < k->ssize; ++i) {
    s_cone_sz = get_sd_cone_size(k->s[i]);
    b[count + i] = s_cone_sz;
  }
  count += k->ssize; /* add ssize here not ssize * (ssize + 1) / 2 */
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
  if (k->ed) {
    c += 3 * k->ed;
  }
  if (k->ep) {
    c += 3 * k->ep;
  }
  if (k->psize) {
    c += 3 * k->psize;
  }
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
  return 0;
}

void SCS(finish_cone)(ScsConeWork *c) {
#ifdef USE_LAPACK
  if (c->Xs) {
    scs_free(c->Xs);
  }
  if (c->Z) {
    scs_free(c->Z);
  }
  if (c->e) {
    scs_free(c->e);
  }
  if (c->work) {
    scs_free(c->work);
  }
#endif
  if (c->cone_boundaries) {
    scs_free(c->cone_boundaries);
  }
  if (c->s) {
    scs_free(c->s);
  }
  if (c) {
    scs_free(c);
  }
}

char *SCS(get_cone_header)(const ScsCone *k) {
  char *tmp = (char *)scs_malloc(sizeof(char) * 512);
  scs_int i, soc_vars, sd_vars;
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
  if (k->ep || k->ed) {
    sprintf(tmp + strlen(tmp), "\t  e: exp vars: %li, dual exp vars: %li\n",
            (long)(3 * k->ep), (long)(3 * k->ed));
  }
  if (k->psize && k->p) {
    sprintf(tmp + strlen(tmp), "\t  p: primal + dual power vars: %li\n",
            (long)(3 * k->psize));
  }
  return tmp;
}

static scs_float exp_newton_one_d(scs_float rho, scs_float y_hat,
                                  scs_float z_hat, scs_float w) {
  scs_float t_prev, t = MAX(w - z_hat, MAX(-z_hat, 1e-9));
  scs_float f = 1., fp = 1.;
  scs_int i;
  for (i = 0; i < EXP_CONE_MAX_ITERS; ++i) {
    t_prev = t;
    f = t * (t + z_hat) / rho / rho - y_hat / rho + log(t / rho) + 1;
    fp = (2 * t + z_hat) / rho / rho + 1 / t;

    t = t - f / fp;

    if (t <= -z_hat) {
      t = -z_hat;
      break;
    } else if (t <= 0) {
      t = 0;
      break;
    } else if (ABS(t - t_prev) < CONE_TOL) {
      break;
    } else if (SQRTF(f * f / fp) < CONE_TOL) {
      break;
    }
  }
  if (i == EXP_CONE_MAX_ITERS) {
    scs_printf("warning: exp cone newton step hit maximum %i iters\n", (int)i);
    scs_printf("rho=%1.5e; y_hat=%1.5e; z_hat=%1.5e; w=%1.5e; f=%1.5e, "
               "fp=%1.5e, t=%1.5e, t_prev= %1.5e\n",
               rho, y_hat, z_hat, w, f, fp, t, t_prev);
  }
  return t + z_hat;
}

static void exp_solve_for_x_with_rho(const scs_float *v, scs_float *x,
                                     scs_float rho, scs_float w) {
  x[2] = exp_newton_one_d(rho, v[1], v[2], w);
  x[1] = (x[2] - v[2]) * x[2] / rho;
  x[0] = v[0] - rho;
}

static scs_float exp_calc_grad(const scs_float *v, scs_float *x, scs_float rho,
                               scs_float w) {
  exp_solve_for_x_with_rho(v, x, rho, w);
  if (x[1] <= 1e-12) {
    return x[0];
  }
  return x[0] + x[1] * log(x[1] / x[2]);
}

static void exp_get_rho_ub(const scs_float *v, scs_float *x, scs_float *ub,
                           scs_float *lb) {
  *lb = 0;
  *ub = 0.125;
  while (exp_calc_grad(v, x, *ub, v[1]) > 0) {
    *lb = *ub;
    (*ub) *= 2;
  }
}

/* project onto the exponential cone, v has dimension *exactly* 3 */
static scs_int proj_exp_cone(scs_float *v) {
  scs_int i;
  scs_float ub, lb, rho, g, x[3];
  scs_float r = v[0], s = v[1], t = v[2];

  /* v in cl(Kexp) */
  if ((s * exp(r / s) - t <= CONE_THRESH && s > 0) ||
      (r <= 0 && s == 0 && t >= 0)) {
    return 0;
  }

  /* -v in Kexp^* */
  if ((r > 0 && r * exp(s / r) + exp(1) * t <= CONE_THRESH) ||
      (r == 0 && s <= 0 && t <= 0)) {
    memset(v, 0, 3 * sizeof(scs_float));
    return 0;
  }

  /* special case with analytical solution */
  if (r < 0 && s < 0) {
    v[1] = 0.0;
    v[2] = MAX(v[2], 0);
    return 0;
  }

  /* iterative procedure to find projection, bisects on dual variable: */
  exp_get_rho_ub(v, x, &ub, &lb); /* get starting upper and lower bounds */
  for (i = 0; i < EXP_CONE_MAX_ITERS; ++i) {
    rho = (ub + lb) / 2; /* halfway between upper and lower bounds */
    g = exp_calc_grad(v, x, rho, x[1]); /* calculates gradient wrt dual var */
    if (g > 0) {
      lb = rho;
    } else {
      ub = rho;
    }
    if (ub - lb < CONE_TOL) {
      break;
    }
  }
#if VERBOSITY > 10
  scs_printf("exponential cone proj iters %i\n", (int)i);
#endif
  if (i == EXP_CONE_MAX_ITERS) {
    scs_printf("warning: exp cone outer step hit maximum %i iters\n", (int)i);
    scs_printf("r=%1.5e; s=%1.5e; t=%1.5e\n", r, s, t);
  }
  v[0] = x[0];
  v[1] = x[1];
  v[2] = x[2];
  return 0;
}

static scs_int set_up_sd_cone_work_space(ScsConeWork *c, const ScsCone *k) {
  scs_int i;
#ifdef USE_LAPACK
  blas_int n_max = 0;
  blas_int neg_one = -1;
  blas_int info = 0;
  scs_float wkopt = 0.0;
#if VERBOSITY > 0
#define _STR_EXPAND(tok) #tok
#define _STR(tok) _STR_EXPAND(tok)
  scs_printf("BLAS(func) = '%s'\n", _STR(BLAS(func)));
#endif
  /* eigenvector decomp workspace */
  for (i = 0; i < k->ssize; ++i) {
    if (k->s[i] > n_max) {
      n_max = (blas_int)k->s[i];
    }
  }
  c->Xs = (scs_float *)scs_calloc(n_max * n_max, sizeof(scs_float));
  c->Z = (scs_float *)scs_calloc(n_max * n_max, sizeof(scs_float));
  c->e = (scs_float *)scs_calloc(n_max, sizeof(scs_float));

  /* workspace query */
  BLAS(syev)
  ("Vectors", "Lower", &n_max, c->Xs, &n_max, SCS_NULL, &wkopt, &neg_one,
   &info);

  if (info != 0) {
    scs_printf("FATAL: syev failure, info = %li\n", (long)info);
    return -1;
  }
  c->lwork = (blas_int)(wkopt + 1); /* +1 for int casting safety */
  c->work = (scs_float *)scs_calloc(c->lwork, sizeof(scs_float));

  if (!c->Xs || !c->Z || !c->e || !c->work) {
    return -1;
  }
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
  return 0;
#endif
}

/* size of X is get_sd_cone_size(n) */
static scs_int proj_semi_definite_cone(scs_float *X, const scs_int n,
                                       ScsConeWork *c) {
/* project onto the positive semi-definite cone */
#ifdef USE_LAPACK
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
  SCS(scale_array)(X, NAN, n);
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
 *    (up to scalar scaling factor which we can ignore due to conic prooperty)
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
      CONE_THRESH + POWF(xh, a) * POWF(yh, (1 - a)) >= rh) {
    return;
  }

  /* -v in K_a^* */
  if (xh <= 0 && yh <= 0 &&
      CONE_THRESH + POWF(-xh, a) * POWF(-yh, 1 - a) >=
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
    if (ABS(f) < CONE_TOL) {
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
      status = proj_semi_definite_cone(&(x[count]), k->s[i], c);
      if (status < 0) {
        return status;
      }
      count += get_sd_cone_size(k->s[i]);
    }
  }

  if (k->ep) { /* doesn't use r_y */
               /*
                * exponential cone is not self dual, if s \in K
                * then y \in K^* and so if K is the primal cone
                * here we project onto K^*, via Moreau
                * \Pi_C^*(y) = y + \Pi_C(-y)
                */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < k->ep; ++i) {
      proj_exp_cone(&(x[count + 3 * i]));
    }
    count += 3 * k->ep;
  }

  /* dual exponential cone */
  if (k->ed) { /* doesn't use r_y */
    /*
     * exponential cone is not self dual, if s \in K
     * then y \in K^* and so if K is the primal cone
     * here we project onto K^*, via Moreau
     * \Pi_C^*(y) = y + \Pi_C(-y)
     */
    scs_int idx;
    scs_float r, s, t;
    SCS(scale_array)(&(x[count]), -1, 3 * k->ed); /* x = -x; */
#ifdef _OPENMP
#pragma omp parallel for private(r, s, t, idx)
#endif
    for (i = 0; i < k->ed; ++i) {
      idx = count + 3 * i;
      r = x[idx];
      s = x[idx + 1];
      t = x[idx + 2];

      proj_exp_cone(&(x[idx]));

      x[idx] -= r;
      x[idx + 1] -= s;
      x[idx + 2] -= t;
    }
    count += 3 * k->ed;
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
  /* project onto OTHER cones */
  return 0;
}

ScsConeWork *SCS(init_cone)(ScsCone *k, scs_int m) {
  ScsConeWork *c = (ScsConeWork *)scs_calloc(1, sizeof(ScsConeWork));
  c->k = k;
  c->m = m;
  c->scaled_cones = 0;
  set_cone_boundaries(k, c);
  c->s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  if (k->ssize && k->s) {
    if (set_up_sd_cone_work_space(c, k) < 0) {
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
