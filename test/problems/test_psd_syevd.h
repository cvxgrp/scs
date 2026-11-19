#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * The PSD projection eigendecomposes each block with syevd (heevd for a
 * complex block), and decomposes it again with syev (heev) when that reports
 * no convergence (cones.c, EIGENDECOMPOSITION OF THE PSD BLOCKS).
 *
 * test_psd_repeated_eigs and test_complex_psd_repeated_eigs project blocks
 * a*I + b*v*v' with |v_i| = 1: n - 1 exactly repeated eigenvalues, the case
 * where syevr was slow and where Accelerate's dsyevr did not converge (#531).
 * The projection is known in closed form. Blocks of different sizes share
 * one cone, so each LAPACK call gets its own block's work-array sizes. The
 * check is relative to the largest entry of the input: a wrong eigenvector,
 * a stale buffer or a fallback that does not restore the block is off by
 * about that much, while the projection is accurate to about 1e-14 of it.
 *
 * test_complex_psd_multiplet_solve solves max t s.t. M - t*I is PSD, with M
 * complex and its smallest eigenvalue repeated, so heevd runs in a solve.
 *
 * A build with the test-only flag SCS_FORCE_PSD_FALLBACK sends every block
 * through the fallback; the projection tests then also check that it ran.
 *
 * test_psd_lapack_int_overflow checks that a block too large for 32-bit
 * LAPACK integers fails at setup instead of overflowing inside LAPACK.
 */

/* largest error, relative to the largest entry of the input, of a projection
 * that is right: far above the rounding (about n * eps), far below a bug */
#define PSD_EIG_TOL VERIFY_TOL(1e-6, 1e-3)

/* svec of the n x n block a*I + b*11' (lower triangle by columns, the
 * off-diagonal entries scaled by sqrt(2)) in x, and in p the svec of its
 * projection onto the PSD cone. The block has the eigenvalue a, n - 1 times,
 * and a + n*b, once. */
static void psd_rank1_block(scs_int n, scs_float a, scs_float b, scs_float *x,
                            scs_float *p) {
  scs_int i, j, idx = 0;
  scs_float lam = a + n * b;
  /* projection: a*(I - 11'/n) if only a >= 0, (lam/n)*11' if only lam >= 0 */
  scs_float pd = (a > 0 ? a * (1 - (scs_float)1. / n) : 0) +
                 (lam > 0 ? lam / n : 0);
  scs_float po = (a > 0 ? -a / n : 0) + (lam > 0 ? lam / n : 0);
  for (j = 0; j < n; ++j) {
    for (i = j; i < n; ++i, ++idx) {
      x[idx] = i == j ? a + b : SQRTF(2.) * b;
      p[idx] = i == j ? pd : SQRTF(2.) * po;
    }
  }
}

/* The complex counterpart: v_k = i^k, so v_j * conj(v_k) is 1, i, -1 or -i,
 * exactly. Lower triangle by columns as in cones.c (unpack_csd_lower): the
 * real diagonal entry, then the real and imaginary parts of the entries
 * below it, scaled by sqrt(2). */
static void csd_rank1_block(scs_int n, scs_float a, scs_float b, scs_float *x,
                            scs_float *p) {
  static const scs_float re[] = {1., 0., -1., 0.}, im[] = {0., 1., 0., -1.};
  scs_int i, j, idx = 0;
  scs_float lam = a + n * b;
  scs_float pd = (a > 0 ? a * (1 - (scs_float)1. / n) : 0) +
                 (lam > 0 ? lam / n : 0);
  scs_float po = (a > 0 ? -a / n : 0) + (lam > 0 ? lam / n : 0);
  for (j = 0; j < n; ++j) {
    x[idx] = a + b;
    p[idx++] = pd;
    for (i = j + 1; i < n; ++i) {
      /* entry (i, j) of v v' is i^(i - j) */
      scs_int r = (i - j) % 4;
      x[idx] = SQRTF(2.) * b * re[r];
      p[idx++] = SQRTF(2.) * po * re[r];
      x[idx] = SQRTF(2.) * b * im[r];
      p[idx++] = SQRTF(2.) * po * im[r];
    }
  }
}

/* Project x with the cone k (m entries, all PSD) and compare it with p,
 * relative to the largest entry of x. Returns 0 if it agrees. */
static const char *psd_syevd_check(ScsCone *k, scs_float *x,
                                   const scs_float *p, scs_int m,
                                   const char *name) {
  ScsConeWork *cw = SCS(init_cone)(k, m);
  scs_float scale = SCS(norm_inf)(x, m), err;
  scs_int ret;
#ifdef SCS_FORCE_PSD_FALLBACK
  scs_int fallbacks = SCS(psd_fallbacks);
#endif
  mu_assert("psd_syevd_check: cone init failed", cw);
  ret = SCS(proj_dual_cone)(x, cw, SCS_NULL, SCS_NULL);
  SCS(finish_cone)(cw);
  mu_assert("psd_syevd_check: PSD projection failed", ret == 0);
  err = SCS(norm_inf_diff)(x, p, m);
  scs_printf("%s: largest error %.2e, %.2e relative to the input\n", name, err,
             err / scale);
  mu_assert_less("psd_syevd_check: projection off", err, PSD_EIG_TOL * scale);
#ifdef SCS_FORCE_PSD_FALLBACK
  mu_assert("psd_syevd_check: SCS_FORCE_PSD_FALLBACK set, but no fallback",
            SCS(psd_fallbacks) > fallbacks);
#endif
  return 0;
}

static const char *test_psd_repeated_eigs(void) {
  /* sizes above and below LAPACK's crossovers (25 for divide and conquer,
   * 32 for blocked tridiagonal reduction); in the second block the repeated
   * eigenvalue is positive, in the others it is negative */
  scs_int s[] = {100, 7, 100};
  scs_float a[] = {3., 1., -3.}, b[] = {-0.25, -0.5, 0.25};
  scs_int i, m = 0, off = 0;
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  scs_float *x, *p;
  const char *fail;
  for (i = 0; i < 3; ++i)
    m += s[i] * (s[i] + 1) / 2;
  x = (scs_float *)scs_calloc(m, sizeof(scs_float));
  p = (scs_float *)scs_calloc(m, sizeof(scs_float));
  for (i = 0; i < 3; ++i) {
    psd_rank1_block(s[i], a[i], b[i], &x[off], &p[off]);
    off += s[i] * (s[i] + 1) / 2;
  }
  k->ssize = 3;
  k->s = s;
  fail = psd_syevd_check(k, x, p, m, "test_psd_repeated_eigs");
  scs_free(x);
  scs_free(p);
  scs_free(k);
  return fail;
}

static const char *test_complex_psd_repeated_eigs(void) {
  scs_int cs[] = {40, 6, 40};
  scs_float a[] = {2., 1., -2.}, b[] = {-0.25, -0.5, 0.25};
  scs_int i, m = 0, off = 0;
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  scs_float *x, *p;
  const char *fail;
  for (i = 0; i < 3; ++i)
    m += cs[i] * cs[i];
  x = (scs_float *)scs_calloc(m, sizeof(scs_float));
  p = (scs_float *)scs_calloc(m, sizeof(scs_float));
  for (i = 0; i < 3; ++i) {
    csd_rank1_block(cs[i], a[i], b[i], &x[off], &p[off]);
    off += cs[i] * cs[i];
  }
  k->cssize = 3;
  k->cs = cs;
  fail = psd_syevd_check(k, x, p, m, "test_complex_psd_repeated_eigs");
  scs_free(x);
  scs_free(p);
  scs_free(k);
  return fail;
}

/* max t s.t. M - t*I is PSD, i.e. min -t, with M = 2*I + 0.25*v*v' complex of
 * size n = 40 (v as in csd_rank1_block): its smallest eigenvalue, 2, is
 * repeated n - 1 times, so the optimal value is -2. */
static const char *test_complex_psd_multiplet_solve(void) {
  const scs_int nb = 40, m = nb * nb;
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_float *b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *p = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *Ax = (scs_float *)scs_calloc(nb, sizeof(scs_float));
  scs_int *Ai = (scs_int *)scs_calloc(nb, sizeof(scs_int));
  scs_int Ap[] = {0, 0}, cs[] = {0};
  scs_float c[] = {-1.}, opt = -2.;
  scs_int j, exitflag;
  const char *fail;

  /* b = vec(M), and A = vec(I): s = b - A*t = vec(M - t*I) */
  csd_rank1_block(nb, 2., 0.25, b, p);
  for (j = 0; j < nb; ++j) {
    Ax[j] = 1.;
    Ai[j] = j < nb - 1 ? j * (2 * nb - j) : m - 1; /* diagonal entries */
  }
  Ap[1] = nb;
  cs[0] = nb;
  k->cs = cs;
  k->cssize = 1;
  d->m = m;
  d->n = 1;
  d->b = b;
  d->c = c;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = 1;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  scs_set_default_settings(stgs);
  stgs->eps_abs = VERIFY_TOL(1e-7, 1e-4);
  stgs->eps_rel = VERIFY_TOL(1e-7, 1e-4);
  exitflag = scs(d, k, stgs, sol, &info);
  scs_printf("test_complex_psd_multiplet_solve: pobj %.10f (optimum %.1f)\n",
             info.pobj, opt);
  mu_assert("test_complex_psd_multiplet_solve: expected SCS_SOLVED",
            exitflag == SCS_SOLVED);
  mu_assert_less("test_complex_psd_multiplet_solve: wrong optimal value",
                 ABS(info.pobj - opt), SOL_TOL * (1 + ABS(opt)));
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(d);
  scs_free(k);
  scs_free(stgs);
  scs_free(b);
  scs_free(p);
  scs_free(Ax);
  scs_free(Ai);
  return fail;
}

/* A block of size n needs LAPACK work arrays of 1 + 6n + 2n^2 entries (real)
 * or 1 + 5n + 2n^2 (complex), which overflow a 32-bit LAPACK integer from
 * n = 32767: setup must fail then, before LAPACK is called. */
static const char *test_psd_lapack_int_overflow(void) {
  scs_int n[] = {32767};
  ScsCone *k;
  ScsConeWork *cw;
  if (SCS_BLAS_INT_MAX > 2147483647.) {
    scs_printf("skipping: LAPACK integers of more than 32 bits\n");
    return mu_skipped;
  }
  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  /* m only sizes a buffer that a failed setup never uses; the true sizes,
   * n(n+1)/2 and n^2, would allocate 4 to 8 GB for nothing */
  k->ssize = 1;
  k->s = n;
  cw = SCS(init_cone)(k, 1);
  mu_assert("test_psd_lapack_int_overflow: real block set up", !cw);
  k->ssize = 0;
  k->s = SCS_NULL;
  k->cssize = 1;
  k->cs = n;
  cw = SCS(init_cone)(k, 1);
  mu_assert("test_psd_lapack_int_overflow: complex block set up", !cw);
  scs_free(k);
  return 0;
}
