#include "cones.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "scs.h"
#include "util.h"

/* Tests for the PSD diagonal-congruence (rank-one svec-diagonal) metric.
 *
 * 1. enforce_cone_boundaries on a PSD block must return exactly
 *    rank-one weights w_ij = delta_i * delta_j within clamp bounds.
 * 2. proj_dual_cone under such a metric must produce a Moreau split
 *    x = out - R^{-1} p with out in K, p in K*, and <out, p> = 0.
 *    Memberships are verified with metric-free Euclidean projections.
 */

#define PSD_N (7)
#define PSD_M ((PSD_N * (PSD_N + 1)) / 2)

static scs_int psd_svec_idx(scs_int i, scs_int j, scs_int n) {
  return j * n - (j * (j - 1)) / 2 + (i - j);
}

/* deterministic LCG so the test has no platform variance */
static scs_float psd_rand(scs_int *state) {
  *state = (scs_int)((1103515245u * (unsigned long)*state + 12345u) &
                     0x7fffffffu);
  return 2.0 * ((scs_float)*state / (scs_float)0x7fffffff) - 1.0;
}

static const char *test_psd_metric(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsConeWork *cw, *eucl;
  scs_int s_dims[1] = {PSD_N};
  scs_float r_y[PSD_M], prof[PSD_M], delta[PSD_N];
  scs_float x[PSD_M], out[PSD_M], p[PSD_M], chk[PSD_M];
  scs_float dot, nrm, dist;
  scs_int i, j, seed = 12345;

  k->ssize = 1;
  k->s = s_dims;
  cw = SCS(init_cone)(k, PSD_M);
  eucl = SCS(init_cone)(k, PSD_M);
  mu_assert("psd_metric: cone init failed", cw && eucl);

  /* --- 1. rank-one fit --- */
  for (i = 0; i < PSD_M; ++i) {
    prof[i] = exp(1.5 * psd_rand(&seed)); /* arbitrary positive profile */
  }
  SCS(enforce_cone_boundaries)(cw, prof, &SCS(mean), 1);
  for (i = 0; i < PSD_N; ++i) {
    delta[i] = SQRTF(prof[psd_svec_idx(i, i, PSD_N)]);
  }
  for (j = 0; j < PSD_N; ++j) {
    for (i = j; i < PSD_N; ++i) {
      scs_float w = prof[psd_svec_idx(i, j, PSD_N)];
      mu_assert("psd_metric: fit not rank-one",
                ABS(w - delta[i] * delta[j]) <= 1e-12 * (1. + ABS(w)));
      mu_assert("psd_metric: fit outside clamp bounds",
                w >= DIAG_SCALE_MULT_MIN * (1. - 1e-12) &&
                    w <= DIAG_SCALE_MULT_MAX * (1. + 1e-12));
    }
  }

  /* --- 2. Moreau split under the fitted metric --- */
  memcpy(r_y, prof, PSD_M * sizeof(scs_float));
  for (i = 0; i < PSD_M; ++i) {
    x[i] = 2.0 * psd_rand(&seed);
  }
  memcpy(out, x, PSD_M * sizeof(scs_float));
  mu_assert("psd_metric: projection failed",
            SCS(proj_dual_cone)(out, cw, SCS_NULL, r_y) == 0);

  /* p = R (out - x) must lie in K* (= K); out must lie in K. */
  for (i = 0; i < PSD_M; ++i) {
    p[i] = r_y[i] * (out[i] - x[i]);
  }
  memcpy(chk, out, PSD_M * sizeof(scs_float));
  SCS(proj_dual_cone)(chk, eucl, SCS_NULL, SCS_NULL);
  dist = SCS(norm_inf_diff)(chk, out, PSD_M);
  nrm = SCS(norm_inf)(out, PSD_M);
  mu_assert("psd_metric: primal part not in cone", dist <= 1e-9 * (1. + nrm));

  memcpy(chk, p, PSD_M * sizeof(scs_float));
  SCS(proj_dual_cone)(chk, eucl, SCS_NULL, SCS_NULL);
  dist = SCS(norm_inf_diff)(chk, p, PSD_M);
  nrm = SCS(norm_inf)(p, PSD_M);
  mu_assert("psd_metric: dual part not in cone", dist <= 1e-9 * (1. + nrm));

  dot = SCS(dot)(out, p, PSD_M);
  mu_assert("psd_metric: Moreau parts not orthogonal",
            ABS(dot) <= 1e-9 * (1. + SCS(norm_inf)(x, PSD_M)));

  /* uniform metric must reproduce the Euclidean projection exactly */
  for (i = 0; i < PSD_M; ++i) {
    r_y[i] = 0.25;
    out[i] = x[i];
    chk[i] = x[i];
  }
  SCS(proj_dual_cone)(out, cw, SCS_NULL, r_y);
  SCS(proj_dual_cone)(chk, eucl, SCS_NULL, SCS_NULL);
  mu_assert("psd_metric: uniform metric deviates from Euclidean",
            SCS(norm_inf_diff)(chk, out, PSD_M) <= 1e-13);

  SCS(finish_cone)(cw);
  SCS(finish_cone)(eucl);
  scs_free(k);
  return 0;
}
