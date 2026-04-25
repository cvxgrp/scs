/*
 * Public Anderson acceleration diagnostic counters.
 */

#ifndef AA_STATS_H_GUARD
#define AA_STATS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs_types.h"

/**
 * Lifetime diagnostics for an Anderson acceleration workspace.
 *
 * Counters accumulate from aa_init and are not cleared by aa_reset, including
 * internal resets after safeguard rejections. Rejection causes are split so a
 * caller can tell what needs tuning.
 */
typedef struct {
  /** Internal AA iteration counter. */
  scs_int iter;
  /** Number of AA updates accepted by aa_apply before safeguarding. */
  scs_int n_accept;
  /** Number of AA rejections due to LAPACK errors. */
  scs_int n_reject_lapack;
  /** Number of AA rejections due to rank-zero reduced systems. */
  scs_int n_reject_rank0;
  /** Number of AA rejections due to non-finite weights. */
  scs_int n_reject_nonfinite;
  /** Number of AA rejections due to the weight-norm cap. */
  scs_int n_reject_weight_cap;
  /** Number of AA steps rejected by safeguarding. */
  scs_int n_safeguard_reject;
  /** Rank of the most recent AA solve. */
  scs_int last_rank;
  /** Weight norm from the most recent AA solve. NaN if no solve was attempted. */
  scs_float last_aa_norm;
  /** Regularization used in the most recent AA solve. */
  scs_float last_regularization;
} AaStats;

#ifdef __cplusplus
}
#endif
#endif
