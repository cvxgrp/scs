#include "normalize.h"

#include "linalg.h"
#include "scs.h"

/* copied from linsys/scs_matrix.c */
#define MIN_NORMALIZATION_FACTOR (1e-4)
#define MAX_NORMALIZATION_FACTOR (1e4)

/* Given D, E in scaling normalize b, c and compute primal / dual scales.
 *
 * Recall that	the normalization routine is performing:
 *
 * [P  A' c]   with   [E  0  0] on both sides (D, E diagonal)
 * [A  0  b]          [0  D  0]
 * [c' b' 0]          [0  0  s]
 *
 * which results in:
 *
 * [ EPE   EA'D  sEc ]
 * [ DAE    0    sDb ]
 * [ sc'E  sb'D   0  ]
 *
 * `s` is incorporated into dual_scale and primal_scale
 *
 */
void SCS(normalize_b_c)(ScsScaling *scal, scs_float *b, scs_float *c) {
  scs_int i;
  scs_float sigma;

  /* scale c */
  for (i = 0; i < scal->n; ++i) {
    c[i] *= scal->E[i];
  }
  /* scale b */
  for (i = 0; i < scal->m; ++i) {
    b[i] *= scal->D[i];
  }

  /* calculate primal and dual scales */
  sigma = MAX(SCS(norm_inf)(c, scal->n), SCS(norm_inf)(b, scal->m));
  sigma = sigma < MIN_NORMALIZATION_FACTOR ? 1.0 : sigma;
  sigma = sigma > MAX_NORMALIZATION_FACTOR ? MAX_NORMALIZATION_FACTOR : sigma;
  sigma = SAFEDIV_POS(1.0, sigma);

  /* Scale b, c */
  SCS(scale_array)(c, sigma, scal->n);
  SCS(scale_array)(b, sigma, scal->m);

  /* We assume that primal_scale = dual_scale, otherwise need to refactorize */
  scal->primal_scale = sigma;
  scal->dual_scale = sigma;
}

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
