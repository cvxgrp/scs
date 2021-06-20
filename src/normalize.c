#include "normalize.h"

#include "linalg.h"
#include "scs.h"

void SCS(normalize_warm_start)(ScsWork *w) {
  scs_int i;
  scs_float *D = w->scal->D;
  scs_float *E = w->scal->E;
  scs_float *x = w->u;
  scs_float *y = &(w->u[w->n]);
  scs_float *s = &(w->v[w->n]);
  for (i = 0; i < w->n; ++i) {
    x[i] *= (E[i] * w->scal->dual_scale);
  }
  for (i = 0; i < w->m; ++i) {
    y[i] /= (D[i] / w->scal->primal_scale);
  }
  for (i = 0; i < w->m; ++i) {
    s[i] *= (D[i] * w->scal->dual_scale);
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
