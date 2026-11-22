/* eigCG deflation, shared by the indirect solvers; see eigcg.h.
 * DEFLATE_VECTORS and EIGCG_WINDOW are compile-time settings in glbopts.h.
 * The work on vectors of length n goes through e->ops; everything else runs
 * here, on the host. */

#include "eigcg.h"

#include <string.h>

/* Cholesky factorization of the Gram matrix W' G W of a harvest's vectors
 * (k x k by column, G[c * k + r]; only the lower triangle, r >= c, is read),
 * in place, that keeps the vectors that are numerically independent of the
 * ones before them (EIGCG_CHOL_REL_PIVOT, eigcg.h). The columns are taken in
 * order; column c's pivot is its diagonal entry less its part in the span
 * of the columns kept before it, and the column is kept if the pivot exceeds
 * EIGCG_CHOL_REL_PIVOT times the diagonal entry, and dropped otherwise. A
 * dropped column takes no part in the later columns' pivots. The pivot never
 * exceeds the diagonal entry (it is the entry less a sum of squares), so a
 * diagonal entry that is not positive, or anything that is not a number,
 * drops the column too. On return idx[0 .. kept) are the kept columns in
 * order, and the leading kept x kept block of G (leading dimension kept)
 * holds the Cholesky factor of their principal submatrix, as
 * eigcg_chol_solve reads it. When no column is dropped, every operation is
 * that of the plain factorization, in the same order. Returns kept. Not
 * static: test/problems/test_deflation.h tests it directly. */
scs_int SCS(eigcg_gram_factor)(scs_float *G, scs_int k, scs_int *idx) {
  scs_int c, r, u, i, j, kept = 0;
  for (c = 0; c < k; ++c) {
    const scs_float g = G[c * k + c];
    for (r = c; r < k; ++r) {
      scs_float acc = G[c * k + r];
      for (u = 0; u < kept; ++u) {
        acc -= G[idx[u] * k + r] * G[idx[u] * k + c];
      }
      G[c * k + r] = acc;
    }
    if (!(G[c * k + c] > EIGCG_CHOL_REL_PIVOT * g)) {
      continue; /* dropped */
    }
    G[c * k + c] = SQRTF(G[c * k + c]);
    for (r = c + 1; r < k; ++r) {
      G[c * k + r] /= G[c * k + c];
    }
    idx[kept++] = c;
  }
  if (kept < k) {
    /* compact in place, by column: idx increases, so every entry moves to
     * an index no larger than its own, and none is overwritten before it
     * is read */
    for (j = 0; j < kept; ++j) {
      for (i = j; i < kept; ++i) {
        G[j * kept + i] = G[idx[j] * k + idx[i]];
      }
    }
  }
  return kept;
}

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

/* solve L L' x = rhs in place given the factor from SCS(eigcg_gram_factor) */
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

/* eigendecomposition of the leading j x j block of T (lda win):
 * eigenvectors into E (lda j), ascending eigenvalues into w */
static scs_int eigcg_eig_leading(ScsEigcg *e, scs_int j, scs_float *E,
                                 scs_float *w) {
  blas_int bn = (blas_int)j, lda = (blas_int)j, lwork = (blas_int)e->lwork,
           info = 0;
  scs_int a, b;
  for (a = 0; a < j; ++a) {
    for (b = 0; b < j; ++b) {
      E[a * j + b] = e->T[a * e->win + b];
    }
  }
  BLAS(syev)("V", "U", &bn, E, &lda, w, e->work, &lwork, &info);
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
    if (nrm > EIGCG_MGS_RANK_TOL) {
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
static void eigcg_restart(ScsEigcg *e, scs_float tlink) {
  scs_int win = e->win, nev = e->max, q, a, t;
  scs_float *E = e->E, *w = e->ev, *S = e->S, *TS = e->TS, *H = e->H,
            *G = e->Gr, *T = e->T;
  scs_float one = 1., zero = 0.;
  blas_int bm, bn2, bk, ld1, ld2, info = 0, lwork = (blas_int)e->lwork;
  /* lowest nev of T_win */
  if (eigcg_eig_leading(e, win, E, w)) {
    e->dead = 1;
    return;
  }
  for (a = 0; a < nev; ++a) {
    memcpy(&(S[a * win]), &(E[a * win]), win * sizeof(scs_float));
  }
  /* lowest nev of the leading (win-1) block, zero-padded last row */
  if (eigcg_eig_leading(e, win - 1, E, w)) {
    e->dead = 1;
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
    e->dead = 1;
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
  BLAS(syev)("V", "U", &bm, H, &ld2, w, e->work, &lwork, &info);
  if (info) {
    e->dead = 1;
    return;
  }
  /* G = S U;  V(:, 0:q) = V(:, 0:win) G, through VS */
  bm = (blas_int)win;
  bk = (blas_int)q;
  BLAS(gemm)("N", "N", &bm, &bn2, &bk, &one, S, &ld1, H, &ld2, &zero, G,
             &ld1);
  if (e->ops->mul_small(e->w, e->VS, e->V, win, G, q) ||
      e->ops->copy(e->w, e->V, e->VS, q)) {
    e->dead = 1;
    e->failed = 1;
    return;
  }
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
  e->j = q;
  e->jd = q;
}

/* begin a harvest: column 0 is the first preconditioned residual */
void SCS(eigcg_start)(ScsEigcg *e, const scs_float *z, scs_float ztr) {
  scs_float isr;
  if (!e->harvest) {
    return;
  }
  e->j = 0;
  e->jd = 0;
  e->step = 0;
  e->have_prev = 0;
  e->dead = 0;
  e->failed = 0;
  memset(e->T, 0, (size_t)e->win * e->win * sizeof(scs_float));
  if (!(ztr > 0.)) {
    e->dead = 1;
    return;
  }
  isr = 1.0 / SQRTF(ztr);
  if (e->ops->scale_copy(e->w, e->V, z, isr)) {
    e->dead = 1;
    e->failed = 1;
    return;
  }
  e->j = 1;
  e->step = 1;
}

/* One observation per CG iteration, at the point where alpha_i, beta_{i+1},
 * z_{i+1} and ztr_{i+1} are all fresh. Completes the pending diagonal entry
 * (1/alpha_i + beta_i/alpha_{i-1}), links the incoming Lanczos vector with
 * sqrt(beta_{i+1})/alpha_i, and appends it with the (-1)^step sign that the
 * positive-offdiagonal tridiagonal convention requires. */
void SCS(eigcg_step)(ScsEigcg *e, scs_float alpha, scs_float beta_next,
                     const scs_float *z, scs_float ztr) {
  scs_int win = e->win, cur;
  scs_float d, tlink, sgn;
  scs_float *T = e->T;
  if (!e->harvest || e->dead) {
    return;
  }
  if (!(alpha > 0.) || !(beta_next >= 0.) || !(ztr > 0.)) {
    e->dead = 1;
    return;
  }
  cur = e->j - 1;
  d = 1. / alpha;
  if (e->have_prev) {
    d += e->pb / e->pa;
  }
  T[cur * win + cur] = d;
  e->jd = e->j;
  tlink = SQRTF(beta_next) / alpha;
  if (e->j == win) {
    eigcg_restart(e, tlink);
    if (e->dead) {
      return;
    }
  } else {
    T[e->j * win + cur] = tlink;
    T[cur * win + e->j] = tlink;
  }
  sgn = (e->step % 2) ? -1.0 / SQRTF(ztr) : 1.0 / SQRTF(ztr);
  if (e->ops->scale_copy(e->w, e->V + (size_t)e->j * e->n, z, sgn)) {
    e->dead = 1;
    e->failed = 1;
    return;
  }
  e->j++;
  e->step++;
  e->pa = alpha;
  e->pb = beta_next;
  e->have_prev = 1;
}

/* After the deep solve: extract the lowest-kv Ritz vectors, form their true
 * images through the operator (kv extra matvecs, counted by the caller),
 * keep those that pass the filter, factor their Gram matrix, and keep those
 * that the factorization keeps (numerically independent of the ones before
 * them). Returns the number of matvecs spent; sets e->failed if one of the
 * solver's ops failed. */
static scs_int eigcg_extract(ScsEigcg *e) {
  const ScsEigcgOps *ops = e->ops;
  scs_int n = e->n, jd = e->jd, kv, a, probes = 0, kept = 0;
  scs_float *E = e->E, *w = e->ev;
  e->count = 0;
  kv = MIN(e->max, jd);
  if (kv <= 0 || e->failed) {
    return 0;
  }
  if (eigcg_eig_leading(e, jd, E, w)) {
    return 0;
  }
  if (ops->mul_small(e->w, e->W, e->V, jd, E, kv)) {
    e->failed = 1;
    return 0;
  }
  /* Keep only columns whose eigenresidual A w = theta Mp w actually
   * holds: deflation coefficients scale like 1/lambda, so a grossly
   * wrong vector injects error amplified by the smallest eigenvalue.
   * Exact vectors pass at machine precision; the tolerance is loose. */
  for (a = 0; a < kv; ++a) {
    scs_float *wv = e->W + (size_t)a * n;
    scs_float *awv = e->AW + (size_t)kept * n;
    scs_float rr = 0., mm = 0., th = w[a];
    if (kept != a) {
      if (ops->copy(e->w, e->W + (size_t)kept * n, wv, 1)) {
        e->failed = 1;
        return probes;
      }
      wv = e->W + (size_t)kept * n;
    }
    probes++;
    if (ops->probe(e->w, wv, awv, th, &rr, &mm)) {
      e->failed = 1;
      return probes;
    }
    if (mm > 0. && rr <= EIGCG_RITZ_RES_TOL * EIGCG_RITZ_RES_TOL * mm) {
      w[kept] = th;
      kept++;
    }
  }
  kv = kept;
  if (kv <= 0) {
    return probes; /* matvecs were still spent probing */
  }
  if (ops->gram(e->w, e->W, e->AW, kv, e->G)) {
    e->failed = 1;
    return probes;
  }
  kept = SCS(eigcg_gram_factor)(e->G, kv, e->idx);
  /* move the kept vectors and their images to the front, in order (each
   * column moves to a lower index, so source and destination differ) */
  for (a = 0; a < kept; ++a) {
    scs_int src = e->idx[a];
    if (src != a &&
        (ops->copy(e->w, e->W + (size_t)a * n, e->W + (size_t)src * n, 1) ||
         ops->copy(e->w, e->AW + (size_t)a * n, e->AW + (size_t)src * n,
                   1))) {
      e->failed = 1;
      return probes;
    }
  }
  e->count = kept;
  return probes;
}

scs_int SCS(eigcg_end_solve)(ScsEigcg *e) {
  scs_int probes;
  if (!e->harvest) {
    return 0;
  }
  probes = eigcg_extract(e);
  e->jd = 0;
  e->harvest = 0;
  if (e->failed) {
    e->count = 0;
    return -1;
  }
  return probes;
}

/* Deflation: Galerkin-project the residual against the harvested
 * eigenvector approximations (the vectors the harvest kept), solve
 * (W'AW) c = W'r via the stored Cholesky factor, and remove the
 * correction. Because W spans a near-invariant subspace of the
 * preconditioned operator, CG does not reintroduce the deflated
 * components -- a single projection here suffices (init-CG). Costs dot
 * products and axpys only. */
scs_int SCS(eigcg_deflate)(ScsEigcg *e, scs_float *x, scs_float *r) {
  const ScsEigcgOps *ops = e->ops;
  scs_int kd = e->count;
  scs_float nrm0, nrm1;
  if (kd <= 0) {
    return 0;
  }
  nrm0 = ops->norm(e->w, r);
  if (ops->coeffs(e->w, e->W, kd, r, e->c)) {
    return -1;
  }
  eigcg_chol_solve(e->G, e->c, kd);
  if (ops->project(e->w, e->W, e->AW, kd, e->c, x, r, 0)) {
    return -1;
  }
  /* For exact eigenvectors this projector is orthogonal and cannot
   * grow the residual; growth means the Ritz error is being amplified
   * by 1/lambda on this right-hand side -- take it back. */
  nrm1 = ops->norm(e->w, r);
  if (!(nrm1 <= EIGCG_UNDO_GROWTH * nrm0)) {
    if (ops->project(e->w, e->W, e->AW, kd, e->c, x, r, 1)) {
      return -1;
    }
  }
  return 0;
}

void SCS(eigcg_free)(ScsEigcg *e) {
  if (e->ops) {
    e->ops->free(e->w, e->W);
    e->ops->free(e->w, e->AW);
    e->ops->free(e->w, e->V);
    e->ops->free(e->w, e->VS);
  }
  scs_free(e->G);
  scs_free(e->c);
  scs_free(e->T);
  scs_free(e->E);
  scs_free(e->ev);
  scs_free(e->S);
  scs_free(e->TS);
  scs_free(e->H);
  scs_free(e->Gr);
  scs_free(e->work);
  scs_free(e->idx);
  e->idx = SCS_NULL;
  e->W = e->AW = e->V = e->VS = SCS_NULL;
  e->G = e->c = e->T = e->E = e->ev = e->S = e->TS = e->H = e->Gr = e->work =
      SCS_NULL;
  e->max = e->count = e->harvest = 0;
}

scs_int SCS(eigcg_init)(ScsEigcg *e, ScsLinSysWork *w, const ScsEigcgOps *ops,
                        scs_int n) {
  scs_int win, nev2;
  memset(e, 0, sizeof(ScsEigcg));
  e->w = w;
  e->ops = ops;
  e->n = n;
  e->max = MIN(DEFLATE_VECTORS, n);
#ifdef SFLOAT
  /* the Lanczos harvest and the small Cholesky sit below the single-
   * precision noise floor; deflation is double-precision only */
  e->max = 0;
#endif
  if (e->max <= 0) {
    e->max = 0;
    return 0;
  }
  win = EIGCG_WINDOW;
  win = MAX(win, 3 * e->max); /* >= nev free slots per cycle */
  win = MAX(win, 2 * e->max + 2);
  e->win = win;
  nev2 = 2 * e->max;
  {
    scs_float wkopt = 0., dummy = 0.;
    blas_int bn = (blas_int)win, lda = (blas_int)win, lw = -1, info = 0;
    BLAS(syev)("V", "U", &bn, &dummy, &lda, &wkopt, &wkopt, &lw, &info);
    e->lwork = info ? 0 : (scs_int)wkopt;
  }
  if (e->lwork <= 0) {
    scs_printf("ERROR: eigCG deflation: LAPACK's syev workspace query "
               "failed\n");
    SCS(eigcg_free)(e);
    return -1;
  }
  /* everything is allocated now, and held until SCS(eigcg_free) (eigcg.h) */
  e->W = ops->alloc(w, n, e->max);
  e->AW = e->W ? ops->alloc(w, n, e->max) : SCS_NULL;
  e->V = e->AW ? ops->alloc(w, n, win) : SCS_NULL;
  e->VS = e->V ? ops->alloc(w, n, nev2) : SCS_NULL;
  e->G = (scs_float *)scs_calloc((size_t)e->max * e->max, sizeof(scs_float));
  e->c = (scs_float *)scs_calloc(e->max, sizeof(scs_float));
  e->T = (scs_float *)scs_calloc((size_t)win * win, sizeof(scs_float));
  e->E = (scs_float *)scs_calloc((size_t)win * win, sizeof(scs_float));
  e->ev = (scs_float *)scs_calloc(win, sizeof(scs_float));
  e->S = (scs_float *)scs_calloc((size_t)win * nev2, sizeof(scs_float));
  e->TS = (scs_float *)scs_calloc((size_t)win * nev2, sizeof(scs_float));
  e->H = (scs_float *)scs_calloc((size_t)nev2 * nev2, sizeof(scs_float));
  e->Gr = (scs_float *)scs_calloc((size_t)win * nev2, sizeof(scs_float));
  e->work = (scs_float *)scs_calloc(e->lwork, sizeof(scs_float));
  e->idx = (scs_int *)scs_calloc(e->max, sizeof(scs_int));
  if (!e->W || !e->AW || !e->V || !e->VS || !e->G || !e->c || !e->T ||
      !e->E || !e->ev || !e->S || !e->TS || !e->H || !e->Gr || !e->work ||
      !e->idx) {
    scs_printf("ERROR: eigCG deflation: could not allocate %.3g MB, "
               "4 x %li + %li vectors of length %li (DEFLATE_VECTORS and "
               "EIGCG_WINDOW in glbopts.h; DEFLATE_VECTORS 0 turns "
               "deflation off)\n",
               (double)((size_t)(2 * nev2 + win) * (size_t)n *
                        sizeof(scs_float)) /
                   1e6,
               (long)e->max, (long)win, (long)n);
    SCS(eigcg_free)(e);
    return -1;
  }
  return 0;
}

#else /* no LAPACK: eigCG needs syev, deflation is off */

scs_int SCS(eigcg_init)(ScsEigcg *e, ScsLinSysWork *w, const ScsEigcgOps *ops,
                        scs_int n) {
  memset(e, 0, sizeof(ScsEigcg));
  e->w = w;
  e->ops = ops;
  e->n = n;
  return 0;
}

void SCS(eigcg_free)(ScsEigcg *e) {
  e->max = e->count = e->harvest = 0;
}

scs_int SCS(eigcg_deflate)(ScsEigcg *e, scs_float *x, scs_float *r) {
  (void)e;
  (void)x;
  (void)r;
  return 0;
}

void SCS(eigcg_start)(ScsEigcg *e, const scs_float *z, scs_float ztr) {
  (void)e;
  (void)z;
  (void)ztr;
}

void SCS(eigcg_step)(ScsEigcg *e, scs_float alpha, scs_float beta_next,
                     const scs_float *z, scs_float ztr) {
  (void)e;
  (void)alpha;
  (void)beta_next;
  (void)z;
  (void)ztr;
}

scs_int SCS(eigcg_end_solve)(ScsEigcg *e) {
  e->harvest = 0;
  return 0;
}

#endif /* USE_LAPACK */

void SCS(eigcg_reset)(ScsEigcg *e) {
  e->count = 0;
}

void SCS(eigcg_begin_solve)(ScsEigcg *e, scs_int cold) {
  e->harvest = (e->max > 0 && cold && e->count == 0);
}
