#include "scs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Two modes, selected by argv[1]:
 *
 *   mismatch (default): preset the WRONG integer width before scs_init;
 *     initialization must fail fast with the interface-layer diagnostic.
 *
 *   gnu: preset the GNU Fortran calling-convention flag together with the
 *     CORRECT integer width (values 2 = LP64|GNU or 3 = ILP64|GNU, as an
 *     MKL-backed NumPy does when it initializes MKL first); initialization
 *     and a full solve must succeed, proving the runtime guard compares
 *     only the width bit. */

#define MKL_INTERFACE_LP64 0
#define MKL_INTERFACE_ILP64 1
#define MKL_INTERFACE_GNU 2

int MKL_Set_Interface_Layer(int);

int main(int argc, char **argv) {
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
  ScsWork *w;
  int gnu_mode = argc > 1 && strcmp(argv[1], "gnu") == 0;
  int layer;
  int actual;
  scs_int exitflag;

#ifdef BLAS64
  layer = gnu_mode ? (MKL_INTERFACE_ILP64 | MKL_INTERFACE_GNU)
                   : MKL_INTERFACE_LP64;
#else
  layer = gnu_mode ? (MKL_INTERFACE_LP64 | MKL_INTERFACE_GNU)
                   : MKL_INTERFACE_ILP64;
#endif

  actual = MKL_Set_Interface_Layer(layer);
  if (actual != layer) {
    fprintf(stderr,
            "failed to pre-set MKL interface layer to %d, actual layer is %d\n",
            layer, actual);
    return 2;
  }

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

  if (gnu_mode) {
    exitflag = scs(&d, &k, &stgs, &sol, &info);
    if (exitflag != SCS_SOLVED) {
      fprintf(stderr, "gnu-layer solve failed with exitflag %d\n",
              (int)exitflag);
      return 4;
    }
    printf("gnu-layer solve succeeded (layer %d)\n", layer);
    free(sol.x);
    free(sol.y);
    free(sol.s);
    return 0;
  }

  w = scs_init(&d, &k, &stgs);
  if (w != SCS_NULL) {
    fprintf(stderr,
            "expected scs_init to fail after forcing the wrong MKL interface "
            "layer\n");
    scs_finish(w);
    return 3;
  }

  return 0;
}
