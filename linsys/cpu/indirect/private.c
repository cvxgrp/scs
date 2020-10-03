#include <limits.h>
#include "private.h"

#define CG_BEST_TOL 1e-12
#define CG_BASE_TOL 1.

char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsMatrix *P,
                              const ScsSettings *stgs) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str,
          "lin-sys:  sparse-indirect\n\t  nnz(A): %li, nnz(P): "
          "%li, cg_rate: %2.2f\n",
          (long)A->p[A->n], P ? (long)P->p[P->n] : 0l, stgs->cg_rate);
  return str;
}

char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str, "lin-sys: avg cg its: %2.2f\n",
          (scs_float)p->tot_cg_its / (info->iter + 1));
  p->tot_cg_its = 0;
  return str;
}

/* M = inv ( diag ( RHO_X * I + P + scale * A'A ) ) */
static void get_preconditioner(const ScsMatrix *A, const ScsMatrix *P,
                               const ScsSettings *stgs, ScsLinSysWork *p) {
  scs_int i, k;
  scs_float *M = p->M;

#if EXTRA_VERBOSE > 0
  scs_printf("getting pre-conditioner\n");
#endif

  for (i = 0; i < A->n; ++i) { /* cols */
    M[i] = stgs->rho_x + stgs->scale * SCS(norm_sq)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]);
    if (P) {
      for (k = P->p[i]; k < P->p[i + 1]; k++) {
        /* diagonal element only */
        if (P->i[k] == i) { /* row == col */
          M[i] += P->x[k];
          break;
        }
      }
    }
    M[i] = 1. / M[i];
  }

#if EXTRA_VERBOSE > 0
  scs_printf("finished getting pre-conditioner\n");
#endif
}

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
#if EXTRA_VERBOSE > 0
  SCS(timer) transpose_timer;
  scs_printf("transposing A\n");
  SCS(tic)(&transpose_timer);
#endif

  z = (scs_int *)scs_calloc(m, sizeof(scs_int));
  for (i = 0; i < Ap[n]; i++) z[Ai[i]]++; /* row counts */
  SCS(cumsum)(Cp, z, m);                  /* row pointers */

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

#if EXTRA_VERBOSE > 0
  scs_printf("finished transposing A, time: %1.2es\n",
             SCS(tocq)(&transpose_timer) / 1e3);
#endif
}

void SCS(free_lin_sys_work)(ScsLinSysWork *p) {
  if (p) {
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

/*y = (RHO_X * I + P + scale * A'A)x */
static void mat_vec(const ScsMatrix *A, const ScsMatrix *P,
                    const ScsSettings *s, ScsLinSysWork *p, const scs_float *x,
                    scs_float *y) {
  scs_float *z = p->tmp;
  memset(z, 0, A->m * sizeof(scs_float)); /* z = 0 */
  memset(y, 0, A->n * sizeof(scs_float)); /* y = 0 */
  if (P) {
    SCS(accum_by_p)(P, p, x, y); /* y = Px */
  }
  SCS(accum_by_a)(A, p, x, z);                 /* z = Ax */
  SCS(scale_array)(z, s->scale, A->m);         /* z = scale * Ax */
  SCS(accum_by_atrans)(A, p, z, y);            /* y += A'z, y = Px + scale * A'Ax */
  SCS(add_scaled_array)(y, x, A->n, s->rho_x); /* y = rho_x * x + Px + scale * A'A x */
}

void SCS(accum_by_atrans)(const ScsMatrix *A, ScsLinSysWork *p,
                          const scs_float *x, scs_float *y) {
  SCS(_accum_by_atrans)(A->n, A->x, A->i, A->p, x, y);
}

void SCS(accum_by_a)(const ScsMatrix *A, ScsLinSysWork *p, const scs_float *x,
                     scs_float *y) {
  SCS(_accum_by_atrans)(p->At->n, p->At->x, p->At->i, p->At->p, x, y);
}

static void apply_pre_conditioner(scs_float *M, scs_float *z, scs_float *r,
                                  scs_int n, scs_float *ipzr) {
  scs_int i;
  *ipzr = 0;
  for (i = 0; i < n; ++i) {
    z[i] = r[i] * M[i];
    *ipzr += z[i] * r[i];
  }
}

ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A, const ScsMatrix *P,
                                      const ScsSettings *stgs) {
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  p->p = (scs_float *)scs_malloc((A->n) * sizeof(scs_float));
  p->r = (scs_float *)scs_malloc((A->n) * sizeof(scs_float));
  p->Gp = (scs_float *)scs_malloc((A->n) * sizeof(scs_float));
  p->tmp = (scs_float *)scs_malloc((A->m) * sizeof(scs_float));

  /* memory for A transpose */
  p->At = (ScsMatrix *)scs_malloc(sizeof(ScsMatrix));
  p->At->m = A->n;
  p->At->n = A->m;
  p->At->i = (scs_int *)scs_malloc((A->p[A->n]) * sizeof(scs_int));
  p->At->p = (scs_int *)scs_malloc((A->m + 1) * sizeof(scs_int));
  p->At->x = (scs_float *)scs_malloc((A->p[A->n]) * sizeof(scs_float));
  transpose(A, p);

  /* preconditioner memory */
  p->z = (scs_float *)scs_malloc((A->n) * sizeof(scs_float));
  p->M = (scs_float *)scs_malloc((A->n) * sizeof(scs_float));
  get_preconditioner(A, P, stgs, p);

  p->tot_cg_its = 0;
  if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At || !p->At->i || !p->At->p ||
      !p->At->x) {
    SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  return p;
}

/* solves (rho_x * I + P + scale * A'A)x = b, s warm start, solution stored in b */
static scs_int pcg(const ScsMatrix *A, const ScsMatrix *P,
                   const ScsSettings *stgs, ScsLinSysWork *pr,
                   const scs_float *s, scs_float *b, scs_int max_its,
                   scs_float tol) {
  scs_int i, n = A->n;
  scs_float ipzr, ipzr_old, alpha;
  scs_float *p = pr->p;   /* cg direction */
  scs_float *Gp = pr->Gp; /* updated CG direction */
  scs_float *r = pr->r;   /* cg residual */
  scs_float *z = pr->z;   /* for preconditioning */
  scs_float *M = pr->M;   /* inverse diagonal preconditioner */

  if (!s) {
    /* take s = 0 */
    /* r = b */
    memcpy(r, b, n * sizeof(scs_float));
    /* b = 0 */
    memset(b, 0, n * sizeof(scs_float));
  } else {
    /* r = Mat * s */
    mat_vec(A, P, stgs, pr, s, r);
    /* r = Mat * s - b */
    SCS(add_scaled_array)(r, b, n, -1);
    /* r = b - Mat * s */
    SCS(scale_array)(r, -1, n);
    /* b = s */
    memcpy(b, s, n * sizeof(scs_float));
  }

  /* check to see if we need to run CG at all */
  if (SCS(norm)(r, n) < MIN(tol, 1e-18)) {
    return 0;
  }

  apply_pre_conditioner(M, z, r, n, &ipzr);
  memcpy(p, z, n * sizeof(scs_float));

  for (i = 0; i < max_its; ++i) {
    mat_vec(A, P, stgs, pr, p, Gp);
    alpha = ipzr / SCS(dot)(p, Gp, n);
    SCS(add_scaled_array)(b, p, n, alpha);
    SCS(add_scaled_array)(r, Gp, n, -alpha);

#if EXTRA_VERBOSE > 0
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, SCS(norm)(r, n),
               (long)i + 1);
#endif
    if (SCS(norm)(r, n) < tol) {
      return i + 1;
    }
    ipzr_old = ipzr;
    apply_pre_conditioner(M, z, r, n, &ipzr);

    SCS(scale_array)(p, ipzr / ipzr_old, n);
    SCS(add_scaled_array)(p, z, n, 1);
  }
  return i;
}

/* solves Mx = b, for x but stores result in b */
/* s contains warm-start (if available) */
/*
 * [x] = [rho_x I + P       A'   ]^{-1} [rx]
 * [y]   [     A        -I / scale]      [ry]
 *
 * becomes:
 *
 * x = (rho_x I + P + scale * A'A)^{-1} (rx + scale * A'ry)
 * y = scale * (Ax - ry)
 *
*/
scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsMatrix *P,
                           const ScsSettings *stgs, ScsLinSysWork *p,
                           scs_float *b, const scs_float *s, scs_int iter) {
  scs_int cg_its, max_iters;
  scs_float cg_tol;

  if (SCS(norm)(b, A->n) <= 1e-18) {
    memset(b, 0, A->n * sizeof(scs_float));
    return 0;
  }

  if (iter < 0) {
    cg_tol = CG_BEST_TOL;
    max_iters = INT_MAX;
  } else {
    cg_tol = SCS(norm)(b, A->n) * CG_BASE_TOL /
             POWF((scs_float)iter + 1, stgs->cg_rate);
    /* set max_its to 3 * n (though in theory n is enough for any tol) */
    max_iters = 3 * A->n;
  }

  /* b = [rx; ry] */
  SCS(scale_array)(&(b[A->n]), stgs->scale, A->m);  /* b[n:] = scale * ry */
  SCS(accum_by_atrans)(A, p, &(b[A->n]), b); /* b[:n] = rx + scale * A'ry */
  /* solves (rho_x I + P + scale A'A)x = b, s warm start, solution stored in b */
  cg_its = pcg(A, P, stgs, p, s, b, max_iters, MAX(cg_tol, CG_BEST_TOL)); /* b[:n] = x */
  SCS(scale_array)(&(b[A->n]), -1. / stgs->scale, A->m);  /* b[n:] = -ry */
  SCS(accum_by_a)(A, p, b, &(b[A->n])); /* b[n:] = Ax - ry */
  SCS(scale_array)(&(b[A->n]), stgs->scale, A->m); /* b[n:] = scale * (Ax - ry) = y */
  if (iter >= 0) {
    p->tot_cg_its += cg_its;
  }

  return 0;
}

void SCS(normalize)(ScsMatrix *A, ScsMatrix *P, const ScsSettings *stgs,
                    const ScsCone *k, ScsScaling *scal, ScsConeWork * c) {
  SCS(_normalize)(A, P, stgs, k, scal, c);
}

void SCS(un_normalize)(ScsMatrix *A, ScsMatrix *P, const ScsSettings *stgs,
                       const ScsScaling *scal) {
  SCS(_un_normalize)(A, P, stgs, scal);
}
