#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * The per-row multipliers of adaptive_diag_scale must not refactorize the
 * matrix for a change that the cone blocks then average away.
 *
 * A PSD block does not keep one multiplier per row: update_scale fits the
 * rank-one form the projection can absorb. A single row of the block can
 * ask for a step beyond the update threshold while the fitted block barely
 * moves. If the change is measured before the fit, that row triggers a new
 * factorization, and a reset of the acceleration and of the cone cache, as
 * often as every RESCALING_MIN_ITERS iterations for a metric that has not
 * changed.
 *
 * The problem is a small Lyapunov SDP,
 *
 *     minimize trace(P)  subject to  A'P + PA + I <= 0,  P >= I,
 *
 * with A an 8 x 8 matrix that is stable by a small margin (largest real
 * part of an eigenvalue: -0.044), which is what makes it slow. There are 36
 * variables and two PSD blocks of 36 rows each.
 *
 * Measured on macOS arm64 with the direct solver over 300 solves with the
 * data perturbed by 1e-13: with the change measured before the fit, 23 to
 * 67 scale updates and 3,125 to 11,325 iterations, about one update per 150
 * iterations; measured on the multipliers that are installed, 2 or 3 scale
 * updates and 1,100 to 1,600 iterations (4 updates with the indirect
 * solver).
 *
 * The test asks for at most 10 updates, or failing that for at most one per
 * 250 iterations. The second form is for a linear solver whose arithmetic
 * is not reproducible from run to run (Apple Accelerate with more than one
 * thread). In about 5,000 runs of this test there, about one in a thousand
 * took a slow path of 5,000 to 15,000 iterations, on which the updates grow
 * with the iterations: 12 to 24 updates, one per 400 to 600 iterations. A
 * count alone would fail such a run although nothing is wrong with it. The
 * old trigger makes one update per 120 to 170 iterations there.
 *
 * The tolerance is beyond single precision, so the test is skipped there.
 */
#define LYAP_N (8)
#define LYAP_P ((LYAP_N * (LYAP_N + 1)) / 2)

/* position of entry (i, j), i >= j, in the scaled lower triangle */
static scs_int lyap_svec_idx(scs_int i, scs_int j) {
  return j * LYAP_N - (j * (j - 1)) / 2 + (i - j);
}

static const char *test_psd_scale_updates(void) {
  static const scs_float Am[LYAP_N][LYAP_N] = {
      {-1.37, 0.22, -0.06, -2.32, 0.43, -2.13, 0.91, 0.61},
      {0.83, -1.09, 0.30, -0.54, -0.31, 1.51, -0.58, -0.23},
      {-0.72, -0.52, -2.23, 0.26, -0.29, -0.35, -0.62, 0.09},
      {-1.34, 0.05, 1.31, -2.69, -0.05, 2.82, -0.99, -1.62},
      {-0.43, 1.81, 2.04, -1.22, -1.35, 0.00, -0.59, -1.95},
      {1.02, -0.86, 0.30, 0.08, 1.87, -0.75, -1.16, -0.26},
      {1.04, -1.33, -0.39, 0.23, -1.62, -1.24, -3.63, 0.29},
      {-0.56, 0.93, 0.11, 0.57, -0.08, 0.75, -0.10, -2.90}};
  /* trace(P) at the solution */
  const scs_float opt = 235.52015;
  const scs_int n = LYAP_P, m = 2 * LYAP_P;
  const scs_float sqrt2 = SQRTF(2.0);

  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag, i, j, r, c, l, nnz = 0;
  scs_int s_dims[2] = {LYAP_N, LYAP_N};
  const char *fail;

  scs_float *Ax = (scs_float *)scs_calloc(m * n, sizeof(scs_float));
  scs_int *Ai = (scs_int *)scs_calloc(m * n, sizeof(scs_int));
  scs_int *Ap = (scs_int *)scs_calloc(n + 1, sizeof(scs_int));
  scs_float *b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *cvec = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float B[LYAP_N][LYAP_N], M[LYAP_N][LYAP_N], col[2 * LYAP_P];

#ifdef SFLOAT
  scs_free(Ax);
  scs_free(Ai);
  scs_free(Ap);
  scs_free(b);
  scs_free(cvec);
  scs_free(sol);
  scs_free(stgs);
  scs_free(d);
  scs_free(k);
  return mu_skipped;
#endif

  /* One variable per entry (i, j), i >= j, of P, in the order of the
   * scaled lower triangle. Its column is the image of the symmetric unit
   * matrix B with B_ij = B_ji = 1: svec(A'B + BA) in the first block
   * (s_1 = -(A'P + PA) - I) and -svec(B) in the second (s_2 = P - I). */
  for (j = 0; j < LYAP_N; ++j) {
    for (i = j; i < LYAP_N; ++i) {
      for (r = 0; r < LYAP_N; ++r) {
        for (c = 0; c < LYAP_N; ++c) {
          B[r][c] = 0.;
        }
      }
      B[i][j] = 1.;
      B[j][i] = 1.;
      for (r = 0; r < LYAP_N; ++r) {
        for (c = 0; c < LYAP_N; ++c) {
          M[r][c] = 0.;
          for (l = 0; l < LYAP_N; ++l) {
            M[r][c] += Am[l][r] * B[l][c] + B[r][l] * Am[l][c];
          }
        }
      }
      for (c = 0; c < LYAP_N; ++c) {
        for (r = c; r < LYAP_N; ++r) {
          scs_float sc = (r == c) ? 1. : sqrt2;
          col[lyap_svec_idx(r, c)] = sc * M[r][c];
          col[LYAP_P + lyap_svec_idx(r, c)] = -sc * B[r][c];
        }
      }
      for (r = 0; r < m; ++r) {
        if (col[r] != 0.) {
          Ax[nnz] = col[r];
          Ai[nnz] = r;
          nnz++;
        }
      }
      Ap[lyap_svec_idx(i, j) + 1] = nnz;
      cvec[lyap_svec_idx(i, j)] = (i == j) ? 1. : 0.;
    }
  }
  /* b = -(svec(I), svec(I)) */
  for (i = 0; i < LYAP_N; ++i) {
    b[lyap_svec_idx(i, i)] = -1.;
    b[LYAP_P + lyap_svec_idx(i, i)] = -1.;
  }

  d->m = m;
  d->n = n;
  d->b = b;
  d->c = cvec;
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d->A->m = m;
  d->A->n = n;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  k->s = s_dims;
  k->ssize = 2;

  scs_set_default_settings(stgs);
  stgs->adaptive_diag_scale = 1; /* the setting under test */
  stgs->eps_abs = 1e-6;
  stgs->eps_rel = 1e-6;

  exitflag = scs(d, k, stgs, sol, &info);

  scs_printf("test_psd_scale_updates: %li iterations, %li scale updates\n",
             (long)info.iter, (long)info.scale_updates);

  mu_assert("test_psd_scale_updates: expected SCS_SOLVED",
            exitflag == SCS_SOLVED);
  mu_assert_less("test_psd_scale_updates: objective ERROR",
                 ABS(info.pobj - opt), 1e-2);
  mu_assert("test_psd_scale_updates: too many scale updates",
            info.scale_updates <= 10 ||
                info.scale_updates * 250 <= info.iter);
  fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(Ax);
  scs_free(Ai);
  scs_free(Ap);
  scs_free(b);
  scs_free(cvec);
  scs_free(k);
  scs_free(stgs);
  scs_free(d);
  return fail;
}
