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

#define EIGCG_WINDOW (60)

#ifdef USE_LAPACK
#include "scs_blas.h"

#ifdef __cplusplus
extern "C" {
#endif

void BLAS(syev)(const char *jobz, const char *uplo, blas_int *n, scs_float *a,
                blas_int *lda, scs_float *w, scs_float *work, blas_int *lwork,
                blas_int *info);
void BLAS(gemm)(const char *transa, const char *transb, blas_int *m,
                blas_int *n, blas_int *k, scs_float *alpha, scs_float *a,
                blas_int *lda, scs_float *b, blas_int *ldb, scs_float *beta,
                scs_float *c, blas_int *ldc);

#ifdef __cplusplus
}
#endif

/* Cholesky of a k x k column-major matrix, lower triangle, in place.
 * Returns nonzero on failure (also catches NaN via the negated compare). */
static scs_int eigcg_chol(scs_float *G, scs_int k) {
  scs_int c, r, t;
  for (c = 0; c < k; ++c) {
    for (r = c; r < k; ++r) {
      scs_float acc = G[c * k + r];
      for (t = 0; t < c; ++t) {
        acc -= G[t * k + r] * G[t * k + c];
      }
      G[c * k + r] = acc;
    }
    if (!(G[c * k + c] > 1e-300)) {
      return 1;
    }
    G[c * k + c] = SQRTF(G[c * k + c]);
    for (r = c + 1; r < k; ++r) {
      G[c * k + r] /= G[c * k + c];
    }
  }
  return 0;
}

/* solve L L' x = rhs in place given the factor from eigcg_chol */
static void eigcg_chol_solve(const scs_float *G, scs_float *x, scs_int k) {
  scs_int c, t;
  for (c = 0; c < k; ++c) {
    for (t = 0; t < c; ++t) {
      x[c] -= G[t * k + c] * x[t];
    }
    x[c] /= G[c * k + c];
  }
  for (c = k - 1; c >= 0; --c) {
    for (t = c + 1; t < k; ++t) {
      x[c] -= G[c * k + t] * x[t];
    }
    x[c] /= G[c * k + c];
  }
}

/* eigendecomposition of the leading j x j block of eig_T (lda eig_win):
 * eigenvectors into E (lda j), ascending eigenvalues into w */
static scs_int eigcg_eig_leading(ScsLinSysWork *pr, scs_int j, scs_float *E,
                                 scs_float *w) {
  blas_int bn = (blas_int)j, lda = (blas_int)j,
           lwork = (blas_int)pr->eig_lwork, info = 0;
  scs_int a, b;
  for (a = 0; a < j; ++a) {
    for (b = 0; b < j; ++b) {
      E[a * j + b] = pr->eig_T[a * pr->eig_win + b];
    }
  }
  BLAS(syev)("V", "U", &bn, E, &lda, w, pr->eig_work, &lwork, &info);
  return (scs_int)info;
}

/* modified Gram-Schmidt on the columns of S, in place with compaction;
 * returns the numerical rank kept */
static scs_int eigcg_mgs(scs_float *S, scs_int rows, scs_int cols) {
  scs_int c, k, t, pass, kept = 0;
  for (c = 0; c < cols; ++c) {
    scs_float *v = &(S[c * rows]);
    scs_float nrm = 0.;
    for (pass = 0; pass < 2; ++pass) {
      for (k = 0; k < kept; ++k) {
        const scs_float *qk = &(S[k * rows]);
        scs_float d = 0.;
        for (t = 0; t < rows; ++t) {
          d += qk[t] * v[t];
        }
        for (t = 0; t < rows; ++t) {
          v[t] -= d * qk[t];
        }
      }
    }
    for (t = 0; t < rows; ++t) {
      nrm += v[t] * v[t];
    }
    nrm = SQRTF(nrm);
    if (nrm > 1e-8) {
      for (t = 0; t < rows; ++t) {
        v[t] /= nrm;
      }
      if (kept != c) {
        memcpy(&(S[kept * rows]), v, rows * sizeof(scs_float));
      }
      kept++;
    }
  }
  return kept;
}

/* Thick restart: compress the full window to the lowest-nev Ritz vectors of
 * this window AND of the previous one (the eigCG doubling trick, which is
 * what makes the Ritz pairs converge despite restarting), re-diagonalise,
 * and leave the arrowhead border that couples the compressed basis to the
 * next incoming Lanczos vector. The CG recurrence is never touched. */
static void eigcg_restart(ScsLinSysWork *pr, scs_float tlink) {
  scs_int win = pr->eig_win, nev = pr->dfl_max, n = pr->n, q, a, t;
  scs_float *E = pr->eig_E, *w = pr->eig_w, *S = pr->eig_S, *TS = pr->eig_TS,
            *H = pr->eig_H, *G = pr->eig_G, *T = pr->eig_T;
  scs_float one = 1., zero = 0.;
  blas_int bm, bn2, bk, ld1, ld2, ld3, info = 0,
           lwork = (blas_int)pr->eig_lwork;
  /* lowest nev of T_win */
  if (eigcg_eig_leading(pr, win, E, w)) {
    pr->eig_dead = 1;
    return;
  }
  for (a = 0; a < nev; ++a) {
    memcpy(&(S[a * win]), &(E[a * win]), win * sizeof(scs_float));
  }
  /* lowest nev of the leading (win-1) block, zero-padded last row */
  if (eigcg_eig_leading(pr, win - 1, E, w)) {
    pr->eig_dead = 1;
    return;
  }
  for (a = 0; a < nev; ++a) {
    for (t = 0; t < win - 1; ++t) {
      S[(nev + a) * win + t] = E[a * (win - 1) + t];
    }
    S[(nev + a) * win + win - 1] = 0.;
  }
  q = eigcg_mgs(S, win, 2 * nev);
  if (q <= 0) {
    pr->eig_dead = 1;
    return;
  }
  /* H = S' T S */
  bm = (blas_int)win;
  bn2 = (blas_int)q;
  bk = (blas_int)win;
  ld1 = (blas_int)win;
  BLAS(gemm)("N", "N", &bm, &bn2, &bk, &one, T, &ld1, S, &ld1, &zero, TS,
             &ld1);
  bm = (blas_int)q;
  ld2 = (blas_int)q;
  BLAS(gemm)("T", "N", &bm, &bn2, &bk, &one, S, &ld1, TS, &ld1, &zero, H,
             &ld2);
  BLAS(syev)("V", "U", &bm, H, &ld2, w, pr->eig_work, &lwork, &info);
  if (info) {
    pr->eig_dead = 1;
    return;
  }
  /* G = S U;  V(:, 0:q) = V(:, 0:win) G */
  bm = (blas_int)win;
  bk = (blas_int)q;
  BLAS(gemm)("N", "N", &bm, &bn2, &bk, &one, S, &ld1, H, &ld2, &zero, G,
             &ld1);
  bm = (blas_int)n;
  bk = (blas_int)win;
  ld3 = (blas_int)n;
  BLAS(gemm)("N", "N", &bm, &bn2, &bk, &one, pr->eig_V, &ld3, G, &ld1, &zero,
             pr->eig_VS, &ld3);
  memcpy(pr->eig_V, pr->eig_VS, (size_t)q * n * sizeof(scs_float));
  /* T = diag(ritz values) + arrowhead border to the incoming vector */
  memset(T, 0, (size_t)win * win * sizeof(scs_float));
  for (a = 0; a < q; ++a) {
    T[a * win + a] = w[a];
  }
  for (a = 0; a < q; ++a) {
    scs_float sa = tlink * G[a * win + (win - 1)];
    T[q * win + a] = sa;
    T[a * win + q] = sa;
  }
  pr->eig_j = q;
  pr->eig_jd = q;
}

/* the window is transient: live only between eigcg_start and
 * eigcg_finish, i.e. during the one deep solve per metric change, so the
 * persistent footprint of deflation is just W and AW (2 k n) */
static void eigcg_free_window(ScsLinSysWork *pr) {
  scs_free(pr->eig_V);
  scs_free(pr->eig_T);
  scs_free(pr->eig_E);
  scs_free(pr->eig_w);
  scs_free(pr->eig_S);
  scs_free(pr->eig_TS);
  scs_free(pr->eig_H);
  scs_free(pr->eig_G);
  scs_free(pr->eig_VS);
  scs_free(pr->eig_work);
  pr->eig_V = pr->eig_T = pr->eig_E = pr->eig_w = pr->eig_S = pr->eig_TS =
      pr->eig_H = pr->eig_G = pr->eig_VS = pr->eig_work = SCS_NULL;
}

static scs_int eigcg_alloc_window(ScsLinSysWork *pr) {
  scs_int n = pr->n, win = pr->eig_win, nev2 = 2 * pr->dfl_max;
  pr->eig_V = (scs_float *)scs_calloc((size_t)n * win, sizeof(scs_float));
  pr->eig_T = (scs_float *)scs_calloc((size_t)win * win, sizeof(scs_float));
  pr->eig_E = (scs_float *)scs_calloc((size_t)win * win, sizeof(scs_float));
  pr->eig_w = (scs_float *)scs_calloc(win, sizeof(scs_float));
  pr->eig_S = (scs_float *)scs_calloc((size_t)win * nev2, sizeof(scs_float));
  pr->eig_TS = (scs_float *)scs_calloc((size_t)win * nev2, sizeof(scs_float));
  pr->eig_H = (scs_float *)scs_calloc((size_t)nev2 * nev2, sizeof(scs_float));
  pr->eig_G = (scs_float *)scs_calloc((size_t)win * nev2, sizeof(scs_float));
  pr->eig_VS = (scs_float *)scs_calloc((size_t)n * nev2, sizeof(scs_float));
  pr->eig_work =
      (scs_float *)scs_calloc(MAX(pr->eig_lwork, 1), sizeof(scs_float));
  if (!pr->eig_V || !pr->eig_T || !pr->eig_E || !pr->eig_w || !pr->eig_S ||
      !pr->eig_TS || !pr->eig_H || !pr->eig_G || !pr->eig_VS ||
      !pr->eig_work) {
    eigcg_free_window(pr);
    return 1;
  }
  return 0;
}

/* begin a harvest: column 0 is the first preconditioned residual */
static void eigcg_start(ScsLinSysWork *pr, const scs_float *z, scs_float ztr) {
  scs_int t, n = pr->n;
  scs_float isr;
  if (!pr->eig_V && eigcg_alloc_window(pr)) {
    pr->eig_dead = 1;
    return;
  }
  pr->eig_j = 0;
  pr->eig_jd = 0;
  pr->eig_step = 0;
  pr->eig_have_prev = 0;
  pr->eig_dead = 0;
  memset(pr->eig_T, 0, (size_t)pr->eig_win * pr->eig_win * sizeof(scs_float));
  if (!(ztr > 0.)) {
    pr->eig_dead = 1;
    return;
  }
  isr = 1.0 / SQRTF(ztr);
  for (t = 0; t < n; ++t) {
    pr->eig_V[t] = isr * z[t];
  }
  pr->eig_j = 1;
  pr->eig_step = 1;
}

/* One observation per CG iteration, at the point where alpha_i, beta_{i+1},
 * z_{i+1} and ztr_{i+1} are all fresh. Completes the pending diagonal entry
 * (1/alpha_i + beta_i/alpha_{i-1}), links the incoming Lanczos vector with
 * sqrt(beta_{i+1})/alpha_i, and appends it with the (-1)^step sign that the
 * positive-offdiagonal tridiagonal convention requires. */
static void eigcg_accum(ScsLinSysWork *pr, scs_float alpha, scs_float beta_next,
                        const scs_float *z, scs_float ztr) {
  scs_int n = pr->n, win = pr->eig_win, cur, t;
  scs_float d, tlink, sgn;
  scs_float *T = pr->eig_T;
  if (pr->eig_dead) {
    return;
  }
  if (!(alpha > 0.) || !(beta_next >= 0.) || !(ztr > 0.)) {
    pr->eig_dead = 1;
    return;
  }
  cur = pr->eig_j - 1;
  d = 1. / alpha;
  if (pr->eig_have_prev) {
    d += pr->eig_pb / pr->eig_pa;
  }
  T[cur * win + cur] = d;
  pr->eig_jd = pr->eig_j;
  tlink = SQRTF(beta_next) / alpha;
  if (pr->eig_j == win) {
    eigcg_restart(pr, tlink);
    if (pr->eig_dead) {
      return;
    }
  } else {
    T[pr->eig_j * win + cur] = tlink;
    T[cur * win + pr->eig_j] = tlink;
  }
  sgn = (pr->eig_step % 2) ? -1.0 / SQRTF(ztr) : 1.0 / SQRTF(ztr);
  {
    scs_float *col = &(pr->eig_V[pr->eig_j * n]);
    for (t = 0; t < n; ++t) {
      col[t] = sgn * z[t];
    }
  }
  pr->eig_j++;
  pr->eig_step++;
  pr->eig_pa = alpha;
  pr->eig_pb = beta_next;
  pr->eig_have_prev = 1;
}

/* After the deep solve: extract the lowest-kv Ritz vectors, form their true
 * images through the operator (kv extra matvecs, counted by the caller),
 * factor the Gram matrix. Returns the number of matvecs spent. */
static scs_int eigcg_extract(ScsLinSysWork *pr) {
  scs_int n = pr->n, jd = pr->eig_jd, kv, a, b, probes;
  scs_float *E = pr->eig_E, *w = pr->eig_w;
  scs_float one = 1., zero = 0.;
  blas_int bm, bn2, bk, ld1, ld2;
  pr->dfl_count = 0;
  kv = MIN(pr->dfl_max, jd);
  if (kv <= 0 || !pr->eig_V) {
    return 0;
  }
  if (eigcg_eig_leading(pr, jd, E, w)) {
    return 0;
  }
  probes = kv;
  bm = (blas_int)n;
  bn2 = (blas_int)kv;
  bk = (blas_int)jd;
  ld1 = (blas_int)n;
  ld2 = (blas_int)jd;
  BLAS(gemm)("N", "N", &bm, &bn2, &bk, &one, pr->eig_V, &ld1, E, &ld2, &zero,
             pr->dfl_p, &ld1);
  {
    /* Keep only columns whose eigenresidual A w = theta Mp w actually
     * holds: deflation coefficients scale like 1/lambda, so a grossly
     * wrong vector injects error amplified by the smallest eigenvalue.
     * Exact vectors pass at machine precision; the tolerance is loose. */
    scs_int kept = 0, t;
    for (a = 0; a < kv; ++a) {
      scs_float *wv = &(pr->dfl_p[a * n]);
      scs_float *awv = &(pr->dfl_Ap[kept * n]);
      scs_float rr = 0., mm = 0., th = w[a];
      if (kept != a) {
        memcpy(&(pr->dfl_p[kept * n]), wv, n * sizeof(scs_float));
        wv = &(pr->dfl_p[kept * n]);
      }
      mat_vec(pr->A, pr->P, pr, wv, awv);
      for (t = 0; t < n; ++t) {
        scs_float mw = th * wv[t] / pr->M[t];
        rr += (awv[t] - mw) * (awv[t] - mw);
        mm += mw * mw;
      }
      if (mm > 0. && rr <= 0.25 * 0.25 * mm) {
        w[kept] = th;
        kept++;
      }
    }
    kv = kept;
    if (kv <= 0) {
      return probes; /* matvecs were still spent probing */
    }
  }
  for (a = 0; a < kv; ++a) {
    for (b = 0; b < kv; ++b) {
      pr->dfl_G[b * kv + a] =
          SCS(dot)(&(pr->dfl_p[a * n]), &(pr->dfl_Ap[b * n]), n);
    }
  }
  if (eigcg_chol(pr->dfl_G, kv)) {
    return probes;
  }
  pr->dfl_count = kv;
  return probes;
}

static scs_int eigcg_finish(ScsLinSysWork *pr) {
  scs_int probes = eigcg_extract(pr);
  eigcg_free_window(pr);
  pr->eig_jd = 0;
  return probes;
}
#endif /* USE_LAPACK */

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

  /* Deflation: Galerkin-project the residual against the harvested
   * eigenvector approximations, solve (W'AW) c = W'r via the stored
   * Cholesky factor, and remove the correction. Because W spans a
   * near-invariant subspace of the preconditioned operator, CG does not
   * reintroduce the deflated components -- a single projection here
   * suffices (init-CG). Costs dot products and axpys only. */
#ifdef USE_LAPACK
  if (pr->dfl_count > 0) {
    scs_int j, k, kd = pr->dfl_count;
    scs_float nrm0 = CG_NORM(r, n), nrm1;
    for (j = 0; j < kd; ++j) {
      pr->dfl_c[j] = SCS(dot)(&(pr->dfl_p[j * n]), r, n);
    }
    eigcg_chol_solve(pr->dfl_G, pr->dfl_c, kd);
    for (j = 0; j < kd; ++j) {
      const scs_float cj = pr->dfl_c[j];
      const scs_float *wj = &(pr->dfl_p[j * n]);
      const scs_float *awj = &(pr->dfl_Ap[j * n]);
      if (cj == 0.) {
        continue;
      }
      for (k = 0; k < n; ++k) {
        b[k] += cj * wj[k];
        r[k] -= cj * awj[k];
      }
    }
    /* For exact eigenvectors this projector is orthogonal and cannot
     * grow the residual; growth means the Ritz error is being amplified
     * by 1/lambda on this right-hand side -- take it back. */
    nrm1 = CG_NORM(r, n);
    if (!(nrm1 <= 1.5 * nrm0)) {
      for (j = 0; j < kd; ++j) {
        const scs_float cj = pr->dfl_c[j];
        const scs_float *wj = &(pr->dfl_p[j * n]);
        const scs_float *awj = &(pr->dfl_Ap[j * n]);
        for (k = 0; k < n; ++k) {
          b[k] -= cj * wj[k];
          r[k] += cj * awj[k];
        }
      }
    }
  }
#endif

  /* check to see if we need to run CG at all */
  if (CG_NORM(r, n) < MAX(tol, 1e-12)) {
    return 0;
  }

  /* z = M r (M is inverse preconditioner) */
  apply_pre_conditioner(z, r, n, pr);
  /* ztr = z'r */
  ztr = SCS(dot)(z, r, n);
  /* p = z */
  memcpy(p, z, n * sizeof(scs_float));
#ifdef USE_LAPACK
  if (pr->dfl_harvest) {
    eigcg_start(pr, z, ztr);
  }
#endif

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
#ifdef USE_LAPACK
    if (pr->dfl_harvest) {
      eigcg_accum(pr, alpha, beta, z, ztr);
    }
#endif
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
  {
    const char *dm = getenv("SCS_DEFLATE");
    char *wm = getenv("SCS_EIGCG_WIN");
    p->dfl_max = dm ? (scs_int)atoi(dm) : DEFLATE_VECTORS;
    p->dfl_max = MIN(p->dfl_max, p->n);
    p->dfl_count = 0;
    p->dfl_harvest = 0;
#ifndef USE_LAPACK
    p->dfl_max = 0; /* eigCG needs syev */
#endif
#ifdef SFLOAT
    /* the Lanczos harvest and the small Cholesky sit below the single-
     * precision noise floor; deflation is double-precision only */
    p->dfl_max = 0;
#endif
    if (p->dfl_max > 0) {
      scs_int win = wm ? (scs_int)atoi(wm) : EIGCG_WINDOW;
      win = MAX(win, 3 * p->dfl_max); /* >= nev free slots per cycle */
      win = MAX(win, 2 * p->dfl_max + 2);
      p->eig_win = win;
      /* only W, AW and the small factor persist; the window itself is
       * allocated per harvest and freed with the deep solve */
      p->dfl_p = (scs_float *)scs_calloc(p->n * p->dfl_max, sizeof(scs_float));
      p->dfl_Ap = (scs_float *)scs_calloc(p->n * p->dfl_max, sizeof(scs_float));
      p->dfl_G =
          (scs_float *)scs_calloc(p->dfl_max * p->dfl_max, sizeof(scs_float));
      p->dfl_c = (scs_float *)scs_calloc(p->dfl_max, sizeof(scs_float));
#ifdef USE_LAPACK
      {
        scs_float wkopt = 0., dummy = 0.;
        blas_int bn = (blas_int)win, lda = (blas_int)win, lw = -1, info = 0;
        BLAS(syev)("V", "U", &bn, &dummy, &lda, &wkopt, &wkopt, &lw, &info);
        p->eig_lwork = info ? 0 : (scs_int)wkopt;
      }
#endif
      if (!p->dfl_p || !p->dfl_Ap || !p->dfl_G || !p->dfl_c ||
          p->eig_lwork <= 0) {
        p->dfl_max = 0;
      }
    }
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
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol) {
  scs_int cg_its, max_iters;
  /* the cold, high-accuracy solve is the one worth harvesting from */
  p->dfl_harvest = (p->dfl_max > 0 && !s && p->dfl_count == 0);

  if (tol <= 0.) {
    scs_printf("Warning: tol = %4f <= 0, likely compiled without setting "
               "INDIRECT flag.\n",
               tol);
  }

  /* use p->tmp here, and in mat_vec, can do both since they don't overlap */
  /* b = [rx; ry] */
  if (CG_NORM(b, p->n + p->m) <= 1e-12) {
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
#ifdef USE_LAPACK
  if (p->dfl_harvest) {
    /* extract the Ritz vectors; the kv image matvecs it spends are added
     * to the iteration count so the matvec accounting stays honest */
    cg_its += eigcg_finish(p);
    p->dfl_harvest = 0;
  }
#endif

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

/* no need to update anything in this case */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  p->diag_r = diag_r; /* this isn't needed but do it to be safe */
  set_preconditioner(p);
  /* the operator and preconditioner changed: the harvested eigenvectors
   * live in the old preconditioned geometry, discard them; the next deep
   * solve (triggered by this same metric change) re-harvests */
  p->dfl_count = 0;
  return 0;
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    scs_free(p->dfl_p);
    scs_free(p->dfl_Ap);
    scs_free(p->dfl_G);
    scs_free(p->dfl_c);
    scs_free(p->eig_V);
    scs_free(p->eig_T);
    scs_free(p->eig_E);
    scs_free(p->eig_w);
    scs_free(p->eig_S);
    scs_free(p->eig_TS);
    scs_free(p->eig_H);
    scs_free(p->eig_G);
    scs_free(p->eig_VS);
    scs_free(p->eig_work);
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
    scs_free(p);
  }
}
