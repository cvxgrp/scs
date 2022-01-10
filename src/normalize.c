#include "normalize.h"

#include "linalg.h"
#include "scs.h"

/* needed for normalizing the warm-start */
void SCS(normalize_sol)(ScsScaling *scal, ScsSolution *sol) {
  scs_int i;
  scs_float *D = scal->D;
  scs_float *E = scal->E;
  for (i = 0; i < scal->n; ++i) {
    sol->x[i] /= (E[i] / scal->dual_scale);
  }
  for (i = 0; i < scal->m; ++i) {
    sol->y[i] /= (D[i] / scal->primal_scale);
  }
  for (i = 0; i < scal->m; ++i) {
    sol->s[i] *= (D[i] * scal->dual_scale);
  }
}

void SCS(un_normalize_sol)(ScsScaling *scal, ScsSolution *sol) {
  scs_int i;
  scs_float *D = scal->D;
  scs_float *E = scal->E;
  for (i = 0; i < scal->n; ++i) {
    sol->x[i] *= (E[i] / scal->dual_scale);
  }
  for (i = 0; i < scal->m; ++i) {
    sol->y[i] *= (D[i] / scal->primal_scale);
  }
  for (i = 0; i < scal->m; ++i) {
    sol->s[i] /= (D[i] * scal->dual_scale);
  }
}

void SCS(un_normalize_primal)(ScsScaling *scal, scs_float *r) {
  scs_int i;
  scs_float *D = scal->D;
  for (i = 0; i < scal->m; ++i) {
    r[i] /= (D[i] * scal->dual_scale);
  }
}

void SCS(un_normalize_dual)(ScsScaling *scal, scs_float *r) {
  scs_int i;
  scs_float *E = scal->E;
  for (i = 0; i < scal->n; ++i) {
    r[i] /= (E[i] * scal->primal_scale);
  }
}
