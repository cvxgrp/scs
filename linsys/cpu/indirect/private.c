#include "private.h"
#include "linsys.h"
#include "util.h"
#include <limits.h>

const char *SCS(get_lin_sys_method)() {
  return "sparse-indirect";
}

/*
char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str, "lin-sys: avg cg its: %2.2f\n",
          (scs_float)p->tot_cg_its / (info->iter + 1));
  p->tot_cg_its = 0;
  return str;
}
*/

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
      for (k = P->p[i]; k < P->p[i + 1]; k++) {
        /* diagonal element only */
        if (P->i[k] == i) { /* row == col */
          /* M_ii += P_ii */
          M[i] += P->x[k];
          break;
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

/* vec -> R_y^{-1} vec */
static void scale_by_r_y_inv(scs_float *vec, ScsLinSysWork *p) {
  scs_int i;
  for (i = 0; i < p->m; ++i) {
    vec[i] /= p->diag_r[p->n + i];
  }
}

/* y += R_x * x  */
static void accum_by_r_x(scs_float *y, const scs_float *x, ScsLinSysWork *p) {
  scs_int i;
  for (i = 0; i < p->n; ++i) {
    y[i] += p->diag_r[i] * x[i];
  }
}

/* we use a different accum_by_a here for speed */
static void accum_by_a(ScsLinSysWork *p, const scs_float *x, scs_float *y) {
  SCS(accum_by_atrans)(p->At, x, y);
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

/* no need to update anything in this case */
void SCS(update_lin_sys_diag_r)(ScsLinSysWork *p, const scs_float *diag_r) {
  p->diag_r = diag_r; /* this isn't needed but do it to be safe */
  set_preconditioner(p);
}

ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A, const ScsMatrix *P,
                                      const scs_float *diag_r) {
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
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
  p->At->m = A->n;
  p->At->n = A->m;
  p->At->i = (scs_int *)scs_calloc((A->p[A->n]), sizeof(scs_int));
  p->At->p = (scs_int *)scs_calloc((A->m + 1), sizeof(scs_int));
  p->At->x = (scs_float *)scs_calloc((A->p[A->n]), sizeof(scs_float));
  transpose(A, p);

  /* preconditioner memory */
  p->diag_r = diag_r;
  p->z = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->M = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  set_preconditioner(p);

  p->tot_cg_its = 0;
  if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At || !p->At->i || !p->At->p ||
      !p->At->x) {
    SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  return p;
}

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

  for (i = 0; i < max_its; ++i) {
    /* Gp = Mat * p */
    mat_vec(A, P, pr, p, Gp);
    /* alpha = z'r / p'G p */
    alpha = ztr / SCS(dot)(p, Gp, n);
    /* b += alpha * p */
    SCS(add_scaled_array)(b, p, n, alpha);
    /* r -= alpha * G p */
    SCS(add_scaled_array)(r, Gp, n, -alpha);

#if VERBOSITY > 3
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, CG_NORM(r, n),
               (long)i + 1);
#endif
    if (CG_NORM(r, n) < tol) {
      return i + 1;
    }
    /* z = M r (M is inverse preconditioner) */
    apply_pre_conditioner(z, r, n, pr);
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
 * [x] = [R_x + P     A' ]^{-1} [rx]
 * [y]   [   A      -R_y ]      [ry]
 *
 * becomes:
 *
 * x = (R_x + P + A' R_y^{-1} A)^{-1} (rx + A' R_y^{-1} ry)
 * y = R_y^{-1} (Ax - ry)
 *
 */
scs_int SCS(solve_lin_sys)(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol) {
  scs_int cg_its, max_iters;

  if (tol <= 0.) {
    scs_printf("Warning: tol = %4f <= 0, likely compiled without setting "
               "INDIRECT flag.\n",
               tol);
  }

  if (CG_NORM(b, p->n + p->m) <= 1e-12) {
    memset(b, 0, (p->n + p->m) * sizeof(scs_float));
    return 0;
  }

  /* use p->tmp here, and in mat_vec, can do both since they don't overlap */
  /* b = [rx; ry] */
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
