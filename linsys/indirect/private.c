#include "private.h"

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsSettings *stgs) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str, "sparse-indirect, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)",
          (long)A->p[A->n], stgs->cg_rate);
  return str;
}

char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info) {
  char *str = (char *)scs_malloc(sizeof(char) * 128);
  sprintf(str,
          "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
          (scs_float)p->tot_cg_its / (info->iter + 1),
          p->total_solve_time / (info->iter + 1) / 1e3);
  p->tot_cg_its = 0;
  p->total_solve_time = 0;
  return str;
}

/* M = inv ( diag ( RHO_X * I + A'A ) ) */
static void get_preconditioner(const ScsMatrix *A, const ScsSettings *stgs,
                               ScsLinSysWork *p) {
  scs_int i;
  scs_float *M = p->M;

#if EXTRA_VERBOSE > 0
  scs_printf("getting pre-conditioner\n");
#endif

  for (i = 0; i < A->n; ++i) {
    M[i] = 1 / (stgs->rho_x +
                SCS(norm_sq)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
    /* M[i] = 1; */
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

/*y = (RHO_X * I + A'A)x */
static void mat_vec(const ScsMatrix *A, const ScsSettings *s, ScsLinSysWork *p,
                    const scs_float *x, scs_float *y) {
  scs_float *tmp = p->tmp;
  memset(tmp, 0, A->m * sizeof(scs_float));
  SCS(accum_by_a)(A, p, x, tmp);
  memset(y, 0, A->n * sizeof(scs_float));
  SCS(accum_by_atrans)(A, p, tmp, y);
  SCS(add_scaled_array)(y, x, A->n, s->rho_x);
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

ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A,
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
  get_preconditioner(A, stgs, p);

  p->total_solve_time = 0;
  p->tot_cg_its = 0;
  if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At || !p->At->i || !p->At->p ||
      !p->At->x) {
    SCS(free_lin_sys_work)(p);
    return SCS_NULL;
  }
  return p;
}

/* solves (I+A'A)x = b, s warm start, solution stored in b */
static scs_int pcg(const ScsMatrix *A, const ScsSettings *stgs,
                   ScsLinSysWork *pr, const scs_float *s, scs_float *b,
                   scs_int max_its, scs_float tol) {
  scs_int i, n = A->n;
  scs_float ipzr, ipzr_old, alpha;
  scs_float *p = pr->p;   /* cg direction */
  scs_float *Gp = pr->Gp; /* updated CG direction */
  scs_float *r = pr->r;   /* cg residual */
  scs_float *z = pr->z;   /* for preconditioning */
  scs_float *M = pr->M;   /* inverse diagonal preconditioner */

  if (s == SCS_NULL) {
    memcpy(r, b, n * sizeof(scs_float));
    memset(b, 0, n * sizeof(scs_float));
  } else {
    mat_vec(A, stgs, pr, s, r);
    SCS(add_scaled_array)(r, b, n, -1);
    SCS(scale_array)(r, -1, n);
    memcpy(b, s, n * sizeof(scs_float));
  }

  /* check to see if we need to run CG at all */
  if (SCS(norm)(r, n) < MIN(tol, 1e-18)) {
    return 0;
  }

  apply_pre_conditioner(M, z, r, n, &ipzr);
  memcpy(p, z, n * sizeof(scs_float));

  for (i = 0; i < max_its; ++i) {
    mat_vec(A, stgs, pr, p, Gp);
    alpha = ipzr / SCS(dot)(p, Gp, n);
    SCS(add_scaled_array)(b, p, n, alpha);
    SCS(add_scaled_array)(r, Gp, n, -alpha);

    if (SCS(norm)(r, n) < tol) {
#if EXTRA_VERBOSE > 0
      scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, SCS(norm)(r, n),
                 (long)i + 1);
#endif
      return i + 1;
    }
    ipzr_old = ipzr;
    apply_pre_conditioner(M, z, r, n, &ipzr);

    SCS(scale_array)(p, ipzr / ipzr_old, n);
    SCS(add_scaled_array)(p, z, n, 1);
  }
  return i;
}

scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsSettings *stgs,
                           ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_int iter) {
  scs_int cg_its;
  SCS(timer) linsys_timer;
  scs_float cg_tol =
      SCS(norm)(b, A->n) *
      (iter < 0 ? CG_BEST_TOL
                : CG_MIN_TOL / POWF((scs_float)iter + 1, stgs->cg_rate));

  SCS(tic)(&linsys_timer);
  /* solves Mx = b, for x but stores result in b */
  /* s contains warm-start (if available) */
  SCS(accum_by_atrans)(A, p, &(b[A->n]), b);
  /* solves (I+A'A)x = b, s warm start, solution stored in b */
  cg_its = pcg(A, stgs, p, s, b, A->n, MAX(cg_tol, CG_BEST_TOL));
  SCS(scale_array)(&(b[A->n]), -1, A->m);
  SCS(accum_by_a)(A, p, b, &(b[A->n]));

  if (iter >= 0) {
    p->tot_cg_its += cg_its;
  }

  p->total_solve_time += SCS(tocq)(&linsys_timer);
#if EXTRA_VERBOSE > 0
  scs_printf("linsys solve time: %1.2es\n", SCS(tocq)(&linsys_timer) / 1e3);
#endif
  return 0;
}

void SCS(normalize_a)(ScsMatrix *A, const ScsSettings *stgs, const ScsCone *k,
                      ScsScaling *scal) {
  SCS(_normalize_a)(A, stgs, k, scal);
}

void SCS(un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                         const ScsScaling *scal) {
  SCS(_un_normalize_a)(A, stgs, scal);
}
