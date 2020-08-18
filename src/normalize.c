#include "normalize.h"

#include "linalg.h"
#include "scs.h"

#define MIN_SCALE (1e-4)

/* Typically l2 equilibration works better than l_inf (Ruiz) */
/* Though more experimentation is needed */
/* #define RUIZ 1 */

#ifdef RUIZ
#define _NORM SCS(norm_inf)
#else
#define _NORM SCS(norm)
#endif

void SCS(normalize_b_c)(ScsWork *w) {
  scs_int i;
  scs_float nm, *D = w->scal->D, *E = w->scal->E, *b = w->b, *c = w->c;
  /* scale b */
  for (i = 0; i < w->m; ++i) {
    b[i] /= D[i];
  }
  if (!w->P || w->scal->dual_scale <= 0.) {
    nm = _NORM(b, w->m);
    w->scal->dual_scale = w->scal->norm_a / MAX(nm, MIN_SCALE);
  }
  /* scale c */
  for (i = 0; i < w->n; ++i) {
    c[i] /= E[i];
  }
  if (!w->P || w->scal->primal_scale <= 0.) {
    nm = _NORM(c, w->n);
    w->scal->primal_scale = w->scal->norm_a / MAX(nm, MIN_SCALE);
  }
  SCS(scale_array)(b, w->scal->dual_scale * w->stgs->scale, w->m);
  SCS(scale_array)(c, w->scal->primal_scale * w->stgs->scale, w->n);
}

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
    y[i] *= (D[i] * w->scal->primal_scale);
  }
  for (i = 0; i < w->m; ++i) {
    s[i] /= (D[i] / (w->scal->dual_scale * w->stgs->scale));
  }
}

void SCS(un_normalize_sol)(ScsWork *w, ScsSolution *sol) {
  scs_int i;
  scs_float *D = w->scal->D;
  scs_float *E = w->scal->E;
  for (i = 0; i < w->n; ++i) {
    sol->x[i] /= (E[i] * w->scal->dual_scale);
  }
  for (i = 0; i < w->m; ++i) {
    sol->y[i] /= (D[i] * w->scal->primal_scale);
  }
  for (i = 0; i < w->m; ++i) {
    sol->s[i] *= D[i] / (w->scal->dual_scale * w->stgs->scale);
  }
}
