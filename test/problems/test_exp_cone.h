#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "scs.h"

/* Forward declare exponential cone projection routine.
 * Implemented in exp_cone.c.
 */
scs_float SCS(proj_pd_exp_cone)(scs_float *v0, scs_int primal);

static scs_int _run_exp_cone_test(scs_float *v0, scs_float *vp_true,
                                  scs_float *vd_true) {
  scs_int success = 1;
  scs_float vp[3], vd[3];
  const scs_float TOL = 1e-6;

  memcpy(vp, v0, 3 * sizeof(scs_float));
  memcpy(vd, v0, 3 * sizeof(scs_float));

  /* inefficient, but just for testing */
  SCS(proj_pd_exp_cone)(vp, 1);
  SCS(proj_pd_exp_cone)(vd, 0);

  scs_printf("*******************************\n");
  scs_printf("v0: (%f, %f, %f)\n", v0[0], v0[1], v0[2]);
  scs_printf("vp: (%f, %f, %f)\n", vp[0], vp[1], vp[2]);
  scs_printf("vp_true: (%f, %f, %f)\n", vp_true[0], vp_true[1], vp_true[2]);
  scs_printf("vd: (%f, %f, %f)\n", vd[0], vd[1], vd[2]);
  scs_printf("vd_true: (%f, %f, %f)\n", vd_true[0], vd_true[1], vd_true[2]);

  success &= (SCS(norm_diff)(vp, vp_true, 3) <= TOL);
  success &= (SCS(norm_diff)(vd, vd_true, 3) <= TOL);
  /* Moreau decomposition holds only for polar */
  /*
  success &= (SCS(dot)(vp, vd, 3) <= TOL);
  success &= (ABS(v0[0] - vp[0] + vd[0]) <= TOL);
  success &= (ABS(v0[1] - vp[1] + vd[1]) <= TOL);
  success &= (ABS(v0[2] - vp[2] + vd[2]) <= TOL);
  */

  if (!success) {
    scs_printf("Failed.\n");
  }

  return success;
}

static const char *test_exp_cone(void) {
  scs_int success = 1;
  scs_int i;
  /* test points */
  scs_float v0[6][3] = {
      {1, 2, 3},
      {0.14814832, 1.04294573, 0.67905585},
      {-0.78301134, 1.82790084, -1.05417044},
      {1.3282585, -0.43277314, 1.7468072},
      {0.67905585, 0.14814832, 1.04294573},
      {0.50210027, 0.12314491, -1.77568921},
  };
  /* primal projections */
  scs_float vp[6][3] = {
      {0.8899428, 1.94041881, 3.06957226},
      {-0.02001571, 0.8709169, 0.85112944},
      {-1.17415616, 0.9567094, 0.280399},
      {0.53160512, 0.2804836, 1.86652094},
      {0.38322814, 0.27086569, 1.11482228},
      {0., 0., 0.},
  };
  /* dual projections */
  scs_float vd[6][3] = {
      {-0., 2., 3.},
      {-0., 1.04294573, 0.67905585},
      {-0.68541419, 1.85424082, 0.01685653},
      {-0.02277033, -0.12164823, 1.75085347},
      {-0., 0.14814832, 1.04294573},
      {-0., 0.12314491, -0.},
  };

  for (i = 0; i < 6; ++i) {
    success &= _run_exp_cone_test(v0[i], vp[i], vd[i]);
  }
  mu_assert("test_exp_cone: Failure", success);
  return 0;
}
