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
  /* data shape */
  int n = 2;
  int m = 3;

  /* Allocate SCS structs */
  ScsCone *k = (ScsCone *)calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)calloc(1, sizeof(ScsSettings));
  ScsSolution *sol = (ScsSolution *)calloc(1, sizeof(ScsSolution));
  ScsInfo *info = (ScsInfo *)calloc(1, sizeof(ScsInfo));

  /* Fill in structs */
  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;
  d->A = &(ScsMatrix){A_x, A_i, A_p, m, n};
  d->P = &(ScsMatrix){P_x, P_i, P_p, n, n};

  /* Cone */
  k->l = m;

  /* Utility to set some default settings */
  scs_set_default_settings(stgs);

  /* Modify tolerances */
  stgs->eps_abs = 1e-9;
  stgs->eps_rel = 1e-9;

  /* Solve! */
  int exitflag = scs(d, k, stgs, sol, info);

  /* Verify that SCS solved the problem */
  printf("SCS solved successfully: %i\n", exitflag == SCS_SOLVED);

  /* Print iterations taken */
  printf("SCS took %i iters\n", info->iter);

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

  return exitflag;
}
