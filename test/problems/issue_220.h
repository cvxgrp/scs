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
 * SCS converges on it, but slowly: the relaxation is degenerate (the Robinson
 * form sits on the boundary of the SOS cone), which is the regime where a
 * first-order method's rate collapses. That is why it looked like a failure
 * to converge in 2022 at 3.2.0. No scale, rho_x, alpha or acceleration
 * setting removes it.
 *
 * The test runs at eps 1e-4, which the iterates reach steadily, in about
 * 1e3 to 4e3 iterations. It must not be tightened to 1e-5. There the
 * residuals sit on a plateau just above the tolerance (the dual residual
 * moves from 2.7e-5 to 2.5e-5 between iteration 1e4 and 2e5) and the solve
 * ends only when a jump happens to land inside it, so the iteration count
 * is decided by rounding: changing `scale` by one part in 1e13 moved it
 * between 6e3 and more than 2e5, and CI runners with different CPUs (and so
 * different BLAS kernels) took 3e4 to more than 2e5, which is the limit
 * stored in the data file.
 *
 * The price of 1e-4 is a less accurate objective: it lands within 2.6e-3 of
 * the exact value, against 8e-4 at 1e-5, hence the tolerance of 5e-3 below.
 * The test guards the data, not a bug.
 */
static const char *issue_220(void) {
  scs_float OPT = 3825. / 4096.; /* exact */
  return _test_prob_from_data_eps("test/problems/issue_220", OPT, 1e-4, 5e-3);
}
