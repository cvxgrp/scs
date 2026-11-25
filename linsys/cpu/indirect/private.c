/* ======================== Includes / Types ======================== */

#include "private.h"

/* ======================== Matrix Helpers ======================== */

static void transpose(const ScsMatrix *A, ScsLinSysWork *p) {
  scs_int *Ci = p->At->i;
  scs_int *Cp = p->At->p;
  scs_float *Cx = p->At->x;
  scs_int m = A->m;
  scs_int n = A->n;

  scs_int *Ap = A->p;
  scs_int *Ai = A->i;
  scs_float *Ax = A->x;

  scs_int i, j, q, *z, c1, c2;
#if VERBOSITY > 0
  SCS(timer) transpose_timer;
  scs_printf("transposing A\n");
  SCS(tic)(&transpose_timer);
#endif

  z = (scs_int *)scs_calloc(m, sizeof(scs_int));
  for (i = 0; i < Ap[n]; i++)
    z[Ai[i]]++;          /* row counts */
  SCS(cumsum)(Cp, z, m); /* row pointers */

  for (j = 0; j < n; j++) {
    c1 = Ap[j];
    c2 = Ap[j + 1];
    for (i = c1; i < c2; i++) {
      q = z[Ai[i]];
      Ci[q] = j; /* place A(i,j) as entry C(j,i) */
      Cx[q] = Ax[i];
      z[Ai[i]]++;
    }
  }
  scs_free(z);

#if VERBOSITY > 0
  scs_printf("finished transposing A, time: %1.2es\n",
             SCS(tocq)(&transpose_timer) / 1e3);
#endif
}

/* Not possible to do this on the fly due to M_ii += a_i' (R_y)^-1 a_i */
/* set M = inv ( diag ( R_x + P + A' R_y^{-1} A ) ) */
static void set_preconditioner(ScsLinSysWork *p) {
  scs_int i, k;
  scs_float *M = p->M;
  const ScsMatrix *A = p->A;
  const ScsMatrix *P = p->P;

#if VERBOSITY > 0
  scs_printf("getting pre-conditioner\n");
#endif

  /* M_ii = (R_x)_i + P_ii + a_i' (R_y)^-1 a_i */
  for (i = 0; i < A->n; ++i) { /* cols */
    /* M_ii = (R_x)_i */
    M[i] = p->diag_r[i];
    /* M_ii += a_i' (R_y)^-1 a_i */
    for (k = A->p[i]; k < A->p[i + 1]; ++k) {
      /* A->i[k] is row of entry k with value A->x[k] */
      M[i] += A->x[k] * A->x[k] / p->diag_r[A->n + A->i[k]];
    }
    if (P) {
      for (k = P->p[i]; k < P->p[i + 1]; ++k) {
        if (P->i[k] == i) {
          M[i] += P->x[k];
        }
      }
    }
    /* finally invert for pre-conditioner */
    M[i] = 1. / M[i];
  }
#if VERBOSITY > 0
  scs_printf("finished getting pre-conditioner\n");
#endif
}

/* we use a different accum_by_a here for speed */
static void accum_by_a(ScsLinSysWork *p, const scs_float *x, scs_float *y) {
  SCS(accum_by_atrans)(p->At, x, y);
}

/* y += R_x * x  */
static void accum_by_r_x(scs_float *y, const scs_float *x, ScsLinSysWork *p) {
  scs_int i;
  for (i = 0; i < p->n; ++i) {
    y[i] += p->diag_r[i] * x[i];
  }
}

/* vec -> R_y^{-1} vec */
static void scale_by_r_y_inv(scs_float *vec, ScsLinSysWork *p) {
  scs_int i;
  for (i = 0; i < p->m; ++i) {
    vec[i] /= p->diag_r[p->n + i];
  }
}

/* y = (R_x + P + A' R_y^{-1} A) x */
static void mat_vec(const ScsMatrix *A, const ScsMatrix *P, ScsLinSysWork *p,
                    const scs_float *x, scs_float *y) {
  scs_float *z = p->tmp;
  memset(z, 0, A->m * sizeof(scs_float)); /* z = 0 */
  memset(y, 0, A->n * sizeof(scs_float)); /* y = 0 */
  if (P) {
    SCS(accum_by_p)(P, x, y); /* y = Px */
  }
  accum_by_a(p, x, z);           /* z = Ax */
  scale_by_r_y_inv(z, p);        /* z = R_y^{-1} A x */
  SCS(accum_by_atrans)(A, z, y); /* y += A'z, y = Px + A' R_y^{-1} Ax */
  /* y = R_x * x + Px + A' R_y^{-1} A * x */
  accum_by_r_x(y, x, p);
}

static void apply_pre_conditioner(scs_float *z, scs_float *r, scs_int n,
                                  ScsLinSysWork *pr) {
  scs_int i;
  scs_float *M = pr->M;
  for (i = 0; i < n; ++i) {
    z[i] = r[i] * M[i];
  }
}

/* ======================== eigCG deflation ======================== */

/* The n-sized operations of the shared eigCG (linsys/eigcg.h), on the
 * host. */

static scs_float *eig_alloc(ScsLinSysWork *p, scs_int rows, scs_int cols) {
  (void)p;
  return (scs_float *)scs_calloc((size_t)rows * cols, sizeof(scs_float));
}

static void eig_free(ScsLinSysWork *p, scs_float *v) {
  (void)p;
  scs_free(v);
}

static scs_int eig_scale_copy(ScsLinSysWork *p, scs_float *dst,
                              const scs_float *src, scs_float a) {
  scs_int t;
  for (t = 0; t < p->n; ++t) {
    dst[t] = a * src[t];
  }
  return 0;
}

static scs_int eig_copy(ScsLinSysWork *p, scs_float *dst, const scs_float *src,
                        scs_int cols) {
  memcpy(dst, src, (size_t)cols * p->n * sizeof(scs_float));
  return 0;
}

#ifdef USE_LAPACK
#include "scs_blas.h"

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(gemm)(const char *transa, const char *transb, blas_int *m,
                blas_int *n, blas_int *k, scs_float *alpha, scs_float *a,
                blas_int *lda, scs_float *b, blas_int *ldb, scs_float *beta,
                scs_float *c, blas_int *ldc);

#ifdef __cplusplus
}
#endif
#endif

static scs_int eig_mul_small(ScsLinSysWork *p, scs_float *C, const scs_float *A,
                             scs_int k, const scs_float *B, scs_int cols) {
#ifdef USE_LAPACK
  scs_float one = 1., zero = 0.;
  blas_int bm = (blas_int)p->n, bn = (blas_int)cols, bk = (blas_int)k,
           ldn = (blas_int)p->n;
  BLAS(gemm)("N", "N", &bm, &bn, &bk, &one, (scs_float *)A, &ldn,
             (scs_float *)B, &bk, &zero, C, &ldn);
  return 0;
#else
  (void)p;
  (void)C;
  (void)A;
  (void)k;
  (void)B;
  (void)cols;
  return -1; /* not reached: deflation is off without LAPACK */
#endif
}

static scs_int eig_probe(ScsLinSysWork *p, const scs_float *v, scs_float *Kv,
                         scs_float theta, scs_float *rr, scs_float *mm) {
  scs_int t;
  scs_float r2 = 0., m2 = 0.;
  mat_vec(p->A, p->P, p, v, Kv);
  for (t = 0; t < p->n; ++t) {
    scs_float mw = theta * v[t] / p->M[t];
    r2 += (Kv[t] - mw) * (Kv[t] - mw);
    m2 += mw * mw;
  }
  *rr = r2;
  *mm = m2;
  return 0;
}

static scs_int eig_gram(ScsLinSysWork *p, const scs_float *W,
                        const scs_float *AW, scs_int k, scs_float *G) {
  scs_int a, b, n = p->n;
  for (a = 0; a < k; ++a) {
    for (b = 0; b < k; ++b) {
      G[b * k + a] = SCS(dot)(W + (size_t)a * n, AW + (size_t)b * n, n);
    }
  }
  return 0;
}

static scs_int eig_coeffs(ScsLinSysWork *p, const scs_float *W, scs_int k,
                          const scs_float *r, scs_float *c) {
  scs_int j, n = p->n;
  for (j = 0; j < k; ++j) {
    c[j] = SCS(dot)(W + (size_t)j * n, r, n);
  }
  return 0;
}

static scs_int eig_project(ScsLinSysWork *p, const scs_float *W,
                           const scs_float *AW, scs_int k, const scs_float *c,
                           scs_float *x, scs_float *r, scs_int undo) {
  scs_int j, t, n = p->n;
  for (j = 0; j < k; ++j) {
    const scs_float cj = c[j];
    const scs_float *wj = W + (size_t)j * n;
    const scs_float *awj = AW + (size_t)j * n;
    if (undo) {
      for (t = 0; t < n; ++t) {
        x[t] -= cj * wj[t];
        r[t] += cj * awj[t];
      }
    } else {
      if (cj == 0.) {
        continue;
      }
      for (t = 0; t < n; ++t) {
        x[t] += cj * wj[t];
        r[t] -= cj * awj[t];
      }
    }
  }
  return 0;
}

static scs_float eig_norm(ScsLinSysWork *p, const scs_float *r) {
  return CG_NORM(r, p->n);
}

static const ScsEigcgOps eig_ops = {eig_alloc,     eig_free,    eig_scale_copy,
                                    eig_copy,      eig_mul_small, eig_probe,
                                    eig_gram,      eig_coeffs,  eig_project,
                                    eig_norm};

/* ======================== Conjugate Gradient Solver ======================== */

/* solves (R_x * I + P + A' R_y^{-1} A)x = b, s warm start, solution in b */
static scs_int pcg(ScsLinSysWork *pr, const scs_float *s, scs_float *b,
                   scs_int max_its, scs_float tol) {
  scs_int i, n = pr->n;
  scs_float ztr, ztr_prev, alpha;
  scs_float *p = pr->p;   /* cg direction */
  scs_float *Gp = pr->Gp; /* updated CG direction */
  scs_float *r = pr->r;   /* cg residual */
  scs_float *z = pr->z;   /* for preconditioning */

  const ScsMatrix *A = pr->A;
  const ScsMatrix *P = pr->P;

  if (!s) {
    /* take s = 0 */
    /* r = b */
    memcpy(r, b, n * sizeof(scs_float));
    /* b = 0 */
    memset(b, 0, n * sizeof(scs_float));
  } else {
    /* r = Mat * s */
    mat_vec(A, P, pr, s, r);
    /* r = Mat * s - b */
    SCS(add_scaled_array)(r, b, n, -1.);
    /* r = b - Mat * s */
    SCS(scale_array)(r, -1., n);
    /* b = s */
    memcpy(b, s, n * sizeof(scs_float));
  }

  /* deflation: project the residual against the harvested eigenvectors
   * (linsys/eigcg.h; the CPU's ops cannot fail) */
  (void)SCS(eigcg_deflate)(&pr->dfl, b, r);

  /* check to see if we need to run CG at all: tol as given (the per-iteration
   * solves' tol is at least CG_BEST_TOL, floored in src/scs.c; the solve for
   * g asks for less) */
  if (CG_NORM(r, n) < tol) {
    return 0;
  }

  /* z = M r (M is inverse preconditioner) */
  apply_pre_conditioner(z, r, n, pr);
  /* ztr = z'r */
  ztr = SCS(dot)(z, r, n);
  /* p = z */
  memcpy(p, z, n * sizeof(scs_float));
  SCS(eigcg_start)(&pr->dfl, z, ztr);

  for (i = 0; i < max_its; ++i) {
    scs_float norm_r, beta;
    scs_int k;
    scs_float *M = pr->M;
    /* Gp = Mat * p */
    mat_vec(A, P, pr, p, Gp);
    /* alpha = z'r / p'G p */
    alpha = ztr / SCS(dot)(p, Gp, n);
    /* b += alpha * p */
    SCS(add_scaled_array)(b, p, n, alpha);
    /* r -= alpha * G p */
    SCS(add_scaled_array)(r, Gp, n, -alpha);

    /* Fuse convergence norm, preconditioner apply, and z'r dot product into
     * one pass over r, saving two separate vector scans per CG iteration. */
    ztr_prev = ztr;
    norm_r = 0.0;
    ztr = 0.0;
    for (k = 0; k < n; ++k) {
      scs_float rk = r[k], zk = rk * M[k], ark = ABS(rk);
      z[k] = zk;
      ztr += zk * rk;
      if (ark > norm_r) norm_r = ark;
    }
#if VERBOSITY > 3
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, norm_r,
               (long)i + 1);
#endif
    if (norm_r < tol) {
      return i + 1;
    }
    if (ztr_prev == 0.) {
      /* preconditioned residual is zero; further CG steps would divide by
       * zero, declare convergence (r must be negligibly small) */
      break;
    }
    /* p = z + beta * p — fuse scale and axpy into one loop */
    beta = ztr / ztr_prev;
    SCS(eigcg_step)(&pr->dfl, alpha, beta, z, ztr);
    for (k = 0; k < n; ++k) {
      p[k] = z[k] + beta * p[k];
    }
  }
  return i;
}

/* ======================== Public API ======================== */

const char *scs_get_lin_sys_method(void) {
  return "sparse-indirect-scs";
}

static scs_int plain_solve(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol);

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p) {
    return SCS_NULL;
  }
  p->A = A;
  p->P = P;
  p->m = A->m;
  p->n = A->n;

  p->p = (scs_float *)scs_calloc((A->n), sizeof(scs_float));
  p->r = (scs_float *)scs_calloc((A->n), sizeof(scs_float));
  p->Gp = (scs_float *)scs_calloc((A->n), sizeof(scs_float));
  p->tmp = (scs_float *)scs_calloc((A->m), sizeof(scs_float));

  /* memory for A transpose */
  p->At = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  p->At->m = A->n;
  p->At->n = A->m;
  p->At->i = (scs_int *)scs_calloc((A->p[A->n]), sizeof(scs_int));
  p->At->p = (scs_int *)scs_calloc((A->m + 1), sizeof(scs_int));
  p->At->x = (scs_float *)scs_calloc((A->p[A->n]), sizeof(scs_float));
  if (!p->At->i || !p->At->p || !p->At->x) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  transpose(A, p);

  /* preconditioner memory */
  p->diag_r = diag_r;
  p->z = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->M = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  if (!p->z || !p->M) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  set_preconditioner(p);

  p->tot_cg_its = 0;
  /* eigCG deflation, with all of its memory (linsys/eigcg.h): setup fails,
   * with a message, if that cannot be allocated */
  if (SCS(eigcg_init)(&p->dfl, p, &eig_ops, p->n) != 0) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  if (SCS(lin_sys_shared_init)(&p->sh, A, P, diag_r, 1, plain_solve,
                               SCS_NULL, SCS_NULL) != 0) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  return p;
}

/* solves Mx = b, for x but stores result in b */
/* s contains warm-start (if available) */
/*
 * [x] = [R_x + P     A' ]^{-1} [rx]
 * [y]   [   A      -R_y ]      [ry]
 *
 * becomes:
 *
 * x = (R_x + P + A' R_y^{-1} A)^{-1} (rx + A' R_y^{-1} ry)
 * y = R_y^{-1} (Ax - ry)
 *
 */
static scs_int plain_solve(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol) {
  scs_int cg_its, max_iters;
  /* the cold solves (s == SCS_NULL: the solve for g, and the refinement
   * steps of scs_solve_lin_sys) are the ones worth harvesting from, the
   * first one only (no set held); warm-started solves do not harvest */
  SCS(eigcg_begin_solve)(&p->dfl, !s);

  if (tol <= 0.) {
    scs_printf("Warning: tol = %4f <= 0, likely compiled without setting "
               "INDIRECT flag.\n",
               tol);
  }

  /* use p->tmp here, and in mat_vec, can do both since they don't overlap */
  /* b = [rx; ry]; x = 0 if b is this small, see CG_ZERO_RHS_TOL */
  if (CG_NORM(b, p->n + p->m) <= MIN(tol, CG_ZERO_RHS_TOL)) {
    memset(b, 0, (p->n + p->m) * sizeof(scs_float));
    return 0;
  }
  /* tmp = ry */
  memcpy(p->tmp, &(b[p->n]), p->m * sizeof(scs_float));
  /* tmp = R_y^{-1} * ry */
  scale_by_r_y_inv(p->tmp, p);
  /* b[:n] = rx + A' R_y^{-1} ry */
  SCS(accum_by_atrans)(p->A, p->tmp, b);
  /* set max_iters to 10 * n (though in theory n is enough for any tol) */
  max_iters = 10 * p->n;
  /* solves (R_x + P + A' R_y^{-1} A)x = b, s warm start, solution stored in
   * b */
  cg_its = pcg(p, s, b, max_iters, tol); /* b[:n] = x */
  /* if it harvested, extract the Ritz vectors; the image matvecs it spends
   * are added to the iteration count so the matvec accounting stays honest
   * (never -1 here: the CPU's ops cannot fail) */
  cg_its += SCS(eigcg_end_solve)(&p->dfl);

  /* b[n:] = -ry */
  SCS(scale_array)(&(b[p->n]), -1., p->m);
  /* b[n:] = Ax - ry */
  accum_by_a(p, b, &(b[p->n]));
  /* b[n:] = R_y^{-1} (Ax - ry) = y */
  scale_by_r_y_inv(&(b[p->n]), p);
  p->tot_cg_its += cg_its;
#if VERBOSITY > 1
  scs_printf("tol %.3e\n", tol);
  scs_printf("cg_its %i\n", (int)cg_its);
#endif
  return 0;
}

/* the warm start, refinement and residual around plain_solve, see linsys.h */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol, scs_int max_refine_steps,
                          scs_float *res) {
  return SCS(lin_sys_shared_solve)(&p->sh, p, b, s, tol, max_refine_steps,
                                   res);
}

/* no need to update anything in this case */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  p->diag_r = diag_r; /* this isn't needed but do it to be safe */
  set_preconditioner(p);
  /* the operator and preconditioner changed: the harvested eigenvectors
   * live in the old preconditioned geometry, discard them; the next deep
   * solve (triggered by this same metric change) re-harvests */
  SCS(eigcg_reset)(&p->dfl);
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    SCS(eigcg_free)(&p->dfl);
    scs_free(p->p);
    scs_free(p->r);
    scs_free(p->Gp);
    scs_free(p->tmp);
    if (p->At) {
      scs_free(p->At->i);
      scs_free(p->At->x);
      scs_free(p->At->p);
      scs_free(p->At);
    }
    scs_free(p->z);
    scs_free(p->M);
    SCS(lin_sys_shared_free)(&p->sh);
    scs_free(p);
  }
}
