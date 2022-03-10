#include "scs.h"    /* SCS API */
#include <stdio.h>  /* printf */
#include <stdlib.h> /* memory allocation */

/* Set up and solve basic qp */
int main(int argc, char **argv) {
  /* Set up the problem data */
  /* A and P must be in compressed sparse column format */
  double P_x[3] = {3., -1., 2.}; /* Upper triangular of P only */
  int P_i[3] = {0, 0, 1};
  int P_p[3] = {0, 1, 3};
  double A_x[4] = {-1., 1., 1., 1.};
  int A_i[4] = {0, 1, 0, 2};
  int A_p[3] = {0, 2, 4};
  double b[3] = {-1., 0.3, -0.5};
  double c[2] = {-1., -1.};
  /* data shapes */
  int n = 2; /* number of variables */
  int m = 3; /* number of constraints */

  /* Allocate SCS structs */
  ScsCone *k = (ScsCone *)calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)calloc(1, sizeof(ScsSolution));
  ScsInfo *info = (ScsInfo *)calloc(1, sizeof(ScsInfo));

  /* Fill in data struct */
  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;
  d->A = &(ScsMatrix){A_x, A_i, A_p, m, n};
  d->P = &(ScsMatrix){P_x, P_i, P_p, n, n};

  /* Cone */
  k->z = 1;
  k->l = 2;

  /* Utility to set default settings */
  scs_set_default_settings(stgs);

  /* Modify tolerances */
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;

  /* Initialize SCS workspace */
  ScsWork *scs_work = scs_init(d, k, stgs);

  /* Solve! */
  int exitflag = scs_solve(scs_work, sol, info, 0);

  /*
   * If we wanted to solve many related problems with different
   * b / c vectors we could update the SCS workspace as follows:
   *
   * int success = scs_update(scs_work, new_b, new_c)
   * int new_exitflag = scs_solve(scs_work, sol, info, 1);
   *
   */

  /* Free SCS workspace */
  scs_finish(scs_work);

  /* Verify that SCS solved the problem */
  printf("SCS solved successfully: %i\n", exitflag == SCS_SOLVED);

  /* Print some info about the solve */
  printf("SCS took %i iters, using the %s linear solver.\n", info->iter,
         info->lin_sys_solver);

  /* Print solution x */
  printf("Optimal solution vector x*:\n");
  for (int i = 0; i < n; ++i) {
    printf("x[%i] = %4f\n", i, sol->x[i]);
  }

  /* Print dual solution y */
  printf("Optimal dual vector y*:\n");
  for (int i = 0; i < m; ++i) {
    printf("y[%i] = %4f\n", i, sol->y[i]);
  }

  /* Free allocated memory */
  free(k);
  free(d);
  free(stgs);
  free(info);
  /* SCS allocates sol->x,y,s if NULL on entry, need to be freed */
  free(sol->x);
  free(sol->y);
  free(sol->s);
  free(sol);
  return 0; /* returning exitflag will set bash exit code to 1 */
}
