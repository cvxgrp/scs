#include "normalize.h"
#include "linalg.h"
#include "scs.h"

#define MIN_SCALE (1e-4)

/* Typically l2 equilibration works better than l1 (Ruiz) */
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
  nm = _NORM(b, w->m);
  w->sc_b = w->scal->norm_a / MAX(nm, MIN_SCALE);
  /* scale c */
  for (i = 0; i < w->n; ++i) {
    c[i] /= E[i];
  }
  if (w->P && w->scal->scale_p > 0) {
    w->sc_c = w->scal->scale_p;
  } else{
    nm = _NORM(c, w->n);
    w->sc_c = w->scal->norm_a / MAX(nm, MIN_SCALE);
  }
  SCS(scale_array)(b, w->sc_b * w->stgs->scale, w->m);
  SCS(scale_array)(c, w->sc_c * w->stgs->scale, w->n);
}

void SCS(normalize_warm_start)(ScsWork *w) {
  scs_int i;
  scs_float *D = w->scal->D;
  scs_float *E = w->scal->E;
  scs_float *x = w->u;
  scs_float *y = &(w->u[w->n]);
  scs_float *s = &(w->v[w->n]);
  for (i = 0; i < w->n; ++i) {
    x[i] *= (E[i] * w->sc_b);
  }
  for (i = 0; i < w->m; ++i) {
    y[i] *= (D[i] * w->sc_c);
  }
  for (i = 0; i < w->m; ++i) {
    s[i] /= (D[i] / (w->sc_b * w->stgs->scale));
  }
}

void SCS(un_normalize_sol)(ScsWork *w, ScsSolution *sol) {
  scs_int i;
  scs_float *D = w->scal->D;
  scs_float *E = w->scal->E;
  for (i = 0; i < w->n; ++i) {
    sol->x[i] /= (E[i] * w->sc_b);
  }
  for (i = 0; i < w->m; ++i) {
    sol->y[i] /= (D[i] * w->sc_c);
  }
  for (i = 0; i < w->m; ++i) {
    sol->s[i] *= D[i] / (w->sc_b * w->stgs->scale);
  }
}
