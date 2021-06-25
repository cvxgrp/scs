#include "normalize.h"

#include "linalg.h"
#include "scs.h"

/* needed for normalizing the warm-start */
void SCS(normalize_sol)(ScsWork *w, ScsSolution *sol) {
  scs_int i;
  scs_float *D = w->scal->D;
  scs_float *E = w->scal->E;
  for (i = 0; i < w->n; ++i) {
    sol->x[i] /= (E[i] / w->scal->dual_scale);
  }
  for (i = 0; i < w->m; ++i) {
    sol->y[i] /= (D[i] / w->scal->primal_scale);
  }
  for (i = 0; i < w->m; ++i) {
    sol->s[i] *= (D[i] * w->scal->dual_scale);
  }
}

void SCS(un_normalize_sol)(ScsWork *w, ScsSolution *sol) {
  scs_int i;
  scs_float *D = w->scal->D;
  scs_float *E = w->scal->E;
  for (i = 0; i < w->n; ++i) {
    sol->x[i] *= (E[i] / w->scal->dual_scale);
  }
  for (i = 0; i < w->m; ++i) {
    sol->y[i] *= (D[i] / w->scal->primal_scale);
  }
  for (i = 0; i < w->m; ++i) {
    sol->s[i] /= (D[i] * w->scal->dual_scale);
  }
}

void SCS(un_normalize_primal)(ScsWork *w, scs_float *r) {
  scs_int i;
  scs_float *D = w->scal->D;
  for (i = 0; i < w->m; ++i) {
    r[i] /= (D[i] * w->scal->dual_scale);
  }
}

void SCS(un_normalize_dual)(ScsWork *w, scs_float *r) {
  scs_int i;
  scs_float *E = w->scal->E;
  for (i = 0; i < w->n; ++i) {
    r[i] /= (E[i] * w->scal->primal_scale);
  }
}
