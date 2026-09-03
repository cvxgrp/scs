/* Minimal consumer of an *installed* SCS package: links scs::scsdir via
 * find_package(scs) and solves a tiny QP. Used by the cmake_consume CI
 * workflow to prove installed packages are usable for both shared and
 * static builds (and that shared installs demand no private build deps). */
#include "scs.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  scs_float P_x[3] = {3., -1., 2.};
  scs_int P_i[3] = {0, 0, 1};
  scs_int P_p[3] = {0, 1, 3};
  scs_float A_x[4] = {-1., 1., 1., 1.};
  scs_int A_i[4] = {0, 1, 0, 2};
  scs_int A_p[3] = {0, 2, 4};
  scs_float b[3] = {-1., 0.3, -0.5};
  scs_float c[2] = {-1., -1.};
  ScsCone k = {0};
  ScsData d = {0};
  ScsSettings stgs = {0};
  ScsSolution sol = {0};
  ScsInfo info;
  scs_int exitflag;

  printf("consuming scs %s\n", scs_version());

  d.m = 3;
  d.n = 2;
  d.b = b;
  d.c = c;
  d.A = &(ScsMatrix){A_x, A_i, A_p, 3, 2};
  d.P = &(ScsMatrix){P_x, P_i, P_p, 2, 2};
  k.z = 1;
  k.l = 2;

  scs_set_default_settings(&stgs);
  stgs.verbose = 1;

  exitflag = scs(&d, &k, &stgs, &sol, &info);
  if (exitflag != SCS_SOLVED) {
    fprintf(stderr, "consume: solve failed with exitflag %d\n", (int)exitflag);
    return 1;
  }
  printf("consume: solved, objective %f\n", info.pobj);
  free(sol.x);
  free(sol.y);
  free(sol.s);
  return 0;
}
