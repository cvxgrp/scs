/*
 * Anderson acceleration (AA) interface for fixed-point iteration.
 * Stores a sliding window of past iterates and computes an accelerated
 * step that can dramatically speed up convergence of ADMM.
 */

#ifndef AA_H_GUARD
#define AA_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "aa_stats.h"
#include "glbopts.h"

typedef scs_float aa_float;
typedef scs_int aa_int;

typedef struct ACCEL_WORK AaWork;

/**
 * Initialize Anderson Acceleration, allocates memory.
 *
 * @param dim               the dimension of the variable for AA
 * @param mem               the memory (number of past iterations used) for AA
 * @param min_len           minimum number of past iterates required before AA
 *                          begins producing updates. Must be >= 1 when
 *                          mem > 0; if min_len exceeds the effective
 *                          memory (min(mem, dim)) it is clamped down,
 *                          mirroring how `mem` is clamped to `dim` for
 *                          rank stability. Set min_len == mem to delay
 *                          AA until the memory is full (stable default
 *                          for large `mem`); set min_len == 1 to start
 *                          extrapolating from the very first residual
 *                          pair (useful when `mem` is large but you
 *                          still want early acceleration). Ignored when
 *                          mem == 0.
 * @param type1             if True use type 1 AA, otherwise use type 2
 * @param regularization    Tikhonov regularization for the AA least-squares
 *                          system. Three modes, selected by sign:
 *                            > 0 : problem-scaled, r = regularization *
 *                                  ||A||_F ||Y||_F. Type-I: 1e-8 works well;
 *                                  Type-II: more stable, 1e-12 often fine.
 *                            < 0 : pinned absolute, r = -regularization
 *                                  (no Frobenius scaling — useful when the
 *                                  problem scale is known).
 *                            = 0 : no regularization.
 *                          Only non-finite values (NaN/Inf) are rejected.
 * @param relaxation        float \in [0,2], mixing parameter (1.0 is vanilla)
 * @param safeguard_factor  factor that controls safeguarding checks
 *                          larger is more aggressive but less stable
 * @param max_weight_norm   float, maximum norm of AA weights
 * @param ir_max_steps      max iterative refinement passes on the γ solve.
 *                          0 disables IR. Each step is O(mem²) and loops
 *                          until the correction stops contracting, so on
 *                          well-conditioned problems only one step runs
 *                          regardless of this cap. Raise it (e.g. 5) for
 *                          ill-conditioned systems where more digits can
 *                          be recovered; lower it for tighter cost bounds.
 * @param verbosity         if greater than 0 prints out various info
 *
 * @return pointer to AA workspace
 *
 */
AaWork *aa_init(aa_int dim, aa_int mem, aa_int min_len, aa_int type1,
                aa_float regularization, aa_float relaxation,
                aa_float safeguard_factor, aa_float max_weight_norm,
                aa_int ir_max_steps, aa_int verbosity);

/**
 * Apply Anderson Acceleration. The usage pattern should be as follows:
 *
 * - for i = 0 .. N:
 *    -  if (i > 0): aa_apply(x, x_prev, a)
 *    -  x_prev = x.copy()
 *    -  x = F(x)
 *    -  aa_safeguard(x, x_prev, a)  (optional but helps stability)
 *
 *  Here F is the map we are trying to find the fixed point for. We put the AA
 *  before the map so that any properties of the map are maintained at the end.
 *  Eg if the map contains a projection onto a set then the output is guaranteed
 *  to be in the set.
 *
 *
 * @param f   output of map at current iteration, overwritten with AA output
 * @param x   input to map at current iteration
 * @param a   workspace from aa_init
 *
 * @return (+ or -) norm of AA weights vector. If positive then update
 *         was accepted and f contains new point, if negative then update was
 *         rejected and f is unchanged
 *
 */
aa_float aa_apply(aa_float *f, const aa_float *x, AaWork *a);

/**
 * Apply safeguarding.
 *
 * This step is optional but can improve stability.
 *
 * @param f_new  output of map after AA step
 * @param x_new  AA output that is input to the map
 * @param a      workspace from aa_init
 *
 * @returns 0 if AA step is accepted otherwise -1, if AA step is rejected then
 *          this overwrites f_new and x_new with previous values
 *
 */
aa_int aa_safeguard(aa_float *f_new, aa_float *x_new, AaWork *a);

/**
 * Finish Anderson Acceleration, clears memory.
 *
 * @param a   AA workspace from aa_init
 *
 */
void aa_finish(AaWork *a);

/**
 * Reset Anderson Acceleration.
 *
 * Resets AA as if at the first iteration, reuses original memory allocations.
 * Does not clear lifetime diagnostic counters; use aa_get_stats after reset
 * to read them, or just re-init the workspace for a clean slate.
 *
 * @param a   AA workspace from aa_init
 *
 */
void aa_reset(AaWork *a);

/**
 * Return lifetime diagnostic counters.
 *
 * Use for post-hoc diagnosis of why AA is or isn't accelerating a
 * given fixed-point iteration — e.g. `n_reject_weight_cap` rising
 * suggests loosening `max_weight_norm` or raising `regularization`;
 * `n_safeguard_reject` rising suggests tuning `safeguard_factor` or
 * `mem`. Do not call with `a == NULL`.
 *
 * @param a    AA workspace from aa_init (must be non-NULL).
 */
AaStats aa_get_stats(const AaWork *a);

#ifdef __cplusplus
}
#endif
#endif
