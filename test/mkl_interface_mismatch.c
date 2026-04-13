#include "scs.h"

#include <stdio.h>
#include <stdlib.h>

#define MKL_INTERFACE_LP64 0
#define MKL_INTERFACE_ILP64 1

int MKL_Set_Interface_Layer(int);

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
  ScsWork *w;
  int wrong_layer;
  int actual;

#ifdef BLAS64
  wrong_layer = MKL_INTERFACE_LP64;
#else
  wrong_layer = MKL_INTERFACE_ILP64;
#endif

  actual = MKL_Set_Interface_Layer(wrong_layer);
  if (actual != wrong_layer) {
    fprintf(stderr,
            "failed to pre-set MKL interface layer to %d, actual layer is %d\n",
            wrong_layer, actual);
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
