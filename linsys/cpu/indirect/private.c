#include <limits.h>
#include "private.h"

#define CG_BASE_TOL (10.)

#ifndef CG_NORM
#define CG_NORM NORM
#endif

char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsMatrix *P,
                              const ScsSettings *stgs) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str,
          "lin-sys:  sparse-indirect\n\t  nnz(A): %li, nnz(P): %li\n",
          (long)A->p[A->n], P ? (long)P->p[P->n] : 0l);
  return str;
}

char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str, "lin-sys: avg cg its: %2.2f\n",
          (scs_float)p->tot_cg_its / (info->iter + 1));
  p->tot_cg_its = 0;
  return str;
}

/* M = inv ( diag ( rho_x * I + P + A' R A ) ) */
static void set_preconditioner(const ScsMatrix *A, const ScsMatrix *P,
                               const ScsSettings *stgs, ScsLinSysWork *p) {
  scs_int i, k;
  scs_float *M = p->M, at_r_a;

#if VERBOSITY > 0
  scs_printf("getting pre-conditioner\n");
#endif

  for (i = 0; i < A->n; ++i) { /* cols */
    M[i] = stgs->rho_x;
    at_r_a = 0.;
    for (k = A->p[i]; k < A->p[i + 1]; ++k) {
      // XXX check that this is correct:
      at_r_a += p->rho_y_vec[A->i[k]] *  A->x[k] * A->x[k];
    }
    M[i] += at_r_a;
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

#if VERBOSITY > 0
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
#if VERBOSITY > 0
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

#if VERBOSITY > 0
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

static void scale_by_diag_r(scs_float *vec, scs_int m, ScsLinSysWork *p) {
  scs_int i;
  for (i = 0; i < m; ++i) {
    vec[i] *= p->rho_y_vec[i];
  }
}

/* y = (rho_x * I + P + A' R A) x */
static void mat_vec(const ScsMatrix *A, const ScsMatrix *P,
                    const ScsSettings *stgs, ScsLinSysWork *p,
                    const scs_float *x, scs_float *y) {
  scs_float *z = p->tmp;
  memset(z, 0, A->m * sizeof(scs_float)); /* z = 0 */
  memset(y, 0, A->n * sizeof(scs_float)); /* y = 0 */
  if (P) {
    SCS(accum_by_p)(P, p, x, y); /* y = Px */
  }
  SCS(accum_by_a)(A, p, x, z); /* z = Ax */
  scale_by_diag_r(z, A->m, p); /* z = R A x */
  SCS(accum_by_atrans)(A, p, z, y); /* y += A'z, y = Px + A' R Ax */
  SCS(add_scaled_array)(y, x, A->n, stgs->rho_x);
  /* y = rho_x * x + Px + A' R A x */
}

void SCS(accum_by_atrans)(const ScsMatrix *A, ScsLinSysWork *p,
                          const scs_float *x, scs_float *y) {
  SCS(_accum_by_atrans)(A->n, A->x, A->i, A->p, x, y);
}

void SCS(accum_by_a)(const ScsMatrix *A, ScsLinSysWork *p, const scs_float *x,
                     scs_float *y) {
  if (p) {
    SCS(_accum_by_atrans)(p->At->n, p->At->x, p->At->i, p->At->p, x, y);
  } else {
    SCS(_accum_by_a)(A->n, A->x, A->i, A->p, x, y, 0);
  }
}

static void apply_pre_conditioner(scs_float *z, scs_float *r,
                                  scs_int n, const ScsSettings *stgs,
                                  ScsLinSysWork *pr) {
  scs_int i;
  scs_float *M = pr->M;
  for (i = 0; i < n; ++i) {
    z[i] = r[i] * M[i];
  }
}

scs_int SCS(should_update_rho_y_vec)(scs_float factor, scs_int iter) {
  return (factor > SQRTF(10.) || factor < 1.0 / SQRTF(10.));
}

/* no need to update anything in this case */
void SCS(update_linsys_rho_y_vec)(const ScsMatrix *A, const ScsMatrix *P,
                                  const ScsSettings *stgs, ScsLinSysWork *p,
                                  scs_float *rho_y_vec) {
  p->rho_y_vec = rho_y_vec;
  set_preconditioner(A, P, stgs, p);
}


ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A, const ScsMatrix *P,
                                      const ScsSettings *stgs,
                                      scs_float *rho_y_vec) {
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
  p->rho_y_vec = rho_y_vec;
  p->z = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->M = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  set_preconditioner(A, P, stgs, p);

  p->tot_cg_its = 0;
  if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At || !p->At->i || !p->At->p ||
      !p->At->x) {
    SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  return p;
}

/* solves (rho_x * I + P + A' R A)x = b, s warm start, solution stored in b */
static scs_int pcg(const ScsMatrix *A, const ScsMatrix *P,
                   const ScsSettings *stgs, ScsLinSysWork *pr,
                   const scs_float *s, scs_float *b, scs_int max_its,
                   scs_float tol) {
  scs_int i, n = A->n;
  scs_float ztr, ztr_prev, alpha;
  scs_float *p = pr->p;   /* cg direction */
  scs_float *Gp = pr->Gp; /* updated CG direction */
  scs_float *r = pr->r;   /* cg residual */
  scs_float *z = pr->z;   /* for preconditioning */

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
    SCS(add_scaled_array)(r, b, n, -1.);
    /* r = b - Mat * s */
    SCS(scale_array)(r, -1., n);
    /* b = s */
    memcpy(b, s, n * sizeof(scs_float));
  }

  /* check to see if we need to run CG at all */
  if (CG_NORM(r, n) < MIN(tol, 1e-12)) {
    return 0;
  }

  /* z = M r (M is inverse preconditioner) */
  apply_pre_conditioner(z, r, n, stgs, pr);
  /* ztr = z'r */
  ztr = SCS(dot)(z, r, n);
  /* p = z */
  memcpy(p, z, n * sizeof(scs_float));

  for (i = 0; i < max_its; ++i) {
    /* Gp = Mat p */
    mat_vec(A, P, stgs, pr, p, Gp);
    /* alpha = z'r / p'G p */
    alpha = ztr / SCS(dot)(p, Gp, n);
    /* b += alpha * p */
    SCS(add_scaled_array)(b, p, n, alpha);
    /* r -= alpha * G p */
    SCS(add_scaled_array)(r, Gp, n, -alpha);

#if VERBOSITY > 1
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, CG_NORM(r, n),
               (long)i + 1);
#endif
    if (CG_NORM(r, n) < tol) {
      return i + 1;
    }
    /* z = M r (M is inverse preconditioner) */
    apply_pre_conditioner(z, r, n, stgs, pr);
    ztr_prev = ztr;
    /* ztr = z'r */
    ztr = SCS(dot)(z, r, n);
    /* p = beta * p, where beta = ztr / ztr_prev */
    SCS(scale_array)(p, ztr / ztr_prev, n);
    /* p = z + beta * p */
    SCS(add_scaled_array)(p, z, n, 1.);
  }
  return i;
}

/* solves Mx = b, for x but stores result in b */
/* s contains warm-start (if available) */
/*
 * [x] = [rho_x I + P     A' ]^{-1} [rx]
 * [y]   [     A        -R^-1]      [ry]
 *
 * R = diag(rho_y_vec)
 *
 * becomes:
 *
 * x = (rho_x I + P + A' R A)^{-1} (rx + A' R ry)
 * y = R (Ax - ry)
 *
 */
scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsMatrix *P,
                           const ScsSettings *stgs, ScsLinSysWork *p,
                           scs_float *b, const scs_float *s, scs_float tol) {
  scs_int cg_its, max_iters;

  if (CG_NORM(b, A->n + A->m) <= 1e-18) {
    memset(b, 0, (A->n + A->m) * sizeof(scs_float));
    return 0;
  }

  /* use p->tmp here, and in mat_vec, can do both since they don't overlap */
  /* b = [rx; ry] */
  /* tmp = ry */
  memcpy(p->tmp, &(b[A->n]), A->m * sizeof(scs_float));
  /* tmp = R * ry */
  scale_by_diag_r(p->tmp, A->m, p);
  /* b[:n] = rx + A' R ry */
  SCS(accum_by_atrans)(A, p, p->tmp, b);
  /* set max_iters to 10 * n (though in theory n is enough for any tol) */
  max_iters = 10 * A->n;
  /* solves (rho_x I + P + A' R A)x = b, s warm start, solution stored in b */
  cg_its = pcg(A, P, stgs, p, s, b, max_iters, tol); /* b[:n] = x */

  /* b[n:] = -ry */
  SCS(scale_array)(&(b[A->n]), -1., A->m);
  /* b[n:] = Ax - ry */
  SCS(accum_by_a)(A, p, b, &(b[A->n]));
  /* b[n:] = R (Ax - ry) = y */
  scale_by_diag_r(&(b[A->n]), A->m, p);
  p->tot_cg_its += cg_its;
#if VERBOSITY > 10
  scs_printf("tol %.3e\n", tol);
  scs_printf("cg_its %i\n", (int)cg_its);
#endif
  return 0;
}

