#include "glbopts.h"
#include "problems/test_prob_from_data_file.h"
#include "scs.h"

/*
 * QP reported in cvxgrp/scs#140: a ~360 variable SOCP that ECOS and Gurobi
 * agree solves to 2.875, but which SCS 2.1.1 could not get near in 100k
 * iterations. The data comes from the `ecos_v_scs.mat` file attached to the
 * issue (github.com/ekrimsk/ECOS_SCS_FILES).
 *
 * The current algorithm solves it in well under a thousand iterations, but
 * only with the diagonal rescaling on: with `adaptive_diag_scale = 0` the
 * primal residual still stalls around 1e-4 for 100k iterations. The data file
 * therefore carries the default settings, and this test guards that default.
 */
static const char *issue_140(void) {
  scs_float OPT = 2.875113; /* from scs at eps 1e-11, matches ECOS to 6 sf */
  return _test_prob_from_data("test/problems/issue_140", OPT);
}
