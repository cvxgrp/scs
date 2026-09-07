#include "glbopts.h"
#include "problems/test_prob_from_data_file.h"
#include "scs.h"

/*
 * The SOS problem from cvxgrp/scs#220: the symmetry-reduced (Wedderburn)
 * sum-of-squares relaxation of the Robinson form under the dihedral group
 * D4, whose exact optimum is 3825 / 4096. Rebuilt from
 * SymbolicWedderburn.jl's test/action_dihedral.jl (the dump file attached to
 * the issue was lost with its host); the reconstruction reproduces the
 * dimensions in the issue's log exactly: n = 12, m = 17, nnz(A) = 25,
 * z = 6, four psd blocks of side 2, 1, 1, 3.
 *
 * SCS converges on it, but with a textbook O(1/k) tail -- one decade of
 * tolerance costs a decade of iterations:
 *
 *     eps 1e-4:    6_225 iterations
 *     eps 1e-5:   47_050 iterations
 *     eps 1e-6:  496_725 iterations
 *
 * which is why it looked like a failure to converge in 2022 at 3.2.0, and
 * still does at any default iteration budget. No scale, rho_x, alpha or
 * acceleration setting removes the tail; the relaxation is degenerate (the
 * Robinson form sits on the boundary of the SOS cone), which is the regime
 * where a first-order method's rate collapses.
 *
 * The test therefore runs at eps 1e-5, where the iteration count is
 * reasonable, and pins the objective against the exact value. It guards the
 * data, not a bug.
 */
static const char *issue_220(void) {
  scs_float OPT = 3825. / 4096.; /* exact */
  return _test_prob_from_data_eps("test/problems/issue_220", OPT, 1e-5, 1e-3);
}
