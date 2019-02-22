#ifndef PUTILS_H_GUARD
#define PUTILS_H_GUARD

#include "amatrix.h"
#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "util.h"

#define PI (3.141592654)
#ifdef DLONG
#ifdef _WIN64
/* this is a Microsoft extension, but also works with min_g_w-w64 */
#define INTRW "%I64d"
#else
#define INTRW "%ld"
#endif
#else
#define INTRW "%i"
#endif

void gen_random_prob_data(scs_int nnz, scs_int col_nnz, ScsData *d, ScsCone *k,
                          ScsSolution *opt_sol);

/* uniform random number in [-1,1] */
static scs_float rand_scs_float(void) {
  return 2 * (((scs_float)rand()) / RAND_MAX) - 1;
}

void gen_random_prob_data(scs_int nnz, scs_int col_nnz, ScsData *d, ScsCone *k,
                          ScsSolution *opt_sol) {
  scs_int n = d->n;
  scs_int m = d->m;
  ScsMatrix *A = d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  scs_float *b = d->b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *c = d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *x = opt_sol->x = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *y = opt_sol->y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *s = opt_sol->s = (scs_float *)scs_calloc(m, sizeof(scs_float));
  /* temporary variables */
  scs_float *z = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_int i, j, r, rn, rm;

  A->i = (scs_int *)scs_calloc(nnz, sizeof(scs_int));
  A->p = (scs_int *)scs_calloc((n + 1), sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(nnz, sizeof(scs_float));
  A->n = d->n;
  A->m = d->m;
  /* y, s >= 0 and y'*s = 0 */
  for (i = 0; i < m; i++) {
    y[i] = z[i] = rand_scs_float();
  }

  SCS(proj_dual_cone)(y, k, SCS_NULL, SCS_NULL, -1);

  for (i = 0; i < m; i++) {
    b[i] = s[i] = y[i] - z[i];
  }

  for (i = 0; i < n; i++) {
    x[i] = rand_scs_float();
  }

  /*
   c = -A'*y
   b = A*x + s
   */
  A->p[0] = 0;
  scs_printf("Generating random matrix:\n");
  for (j = 0; j < n; j++) { /* column */
    if (j * 100 % n == 0 && (j * 100 / n) % 10 == 0) {
      scs_printf("%ld%%\n", (long)(j * 100 / n));
    }
    r = 0;
    for (i = 0; i < m && r < col_nnz; ++i) {
      /* generate a unique sorted array via Knuths alg */
      rn = m - i;
      rm = col_nnz - r;
      if ((rand() % rn) < rm) {
        A->x[r + j * col_nnz] = rand_scs_float();
        A->i[r + j * col_nnz] = i;
        b[i] += A->x[r + j * col_nnz] * x[j];
        c[j] -= A->x[r + j * col_nnz] * y[i];
        r++;
      }
    }
    A->p[j + 1] = (j + 1) * col_nnz;
  }
  scs_printf("done\n");
  scs_free(z);
}

#endif
