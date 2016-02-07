#ifndef PUTILS_H_GUARD
#define PUTILS_H_GUARD

#include "scs.h"
#include "linsys/amatrix.h"

#define PI (3.141592654)
#ifdef DLONG
#ifdef _WIN64
/* this is a Microsoft extension, but also works with MinGW-w64 */
#define INTRW "%I64d"
#else
#define INTRW "%ld"
#endif
#else
#define INTRW "%i"
#endif

#ifndef FLOAT
#define FLOATRW "%lf"
#else
#define FLOATRW "%f"
#endif

/* uniform random number in [-1,1] */
scs_float rand_scs_float(void) {
    return 2 * (((scs_float)rand()) / RAND_MAX) - 1;
}

/* normal random var */
static scs_float U, V;
static scs_int phase = 0;
scs_float rand_gauss(void) {
    scs_float Z;
    if (phase == 0) {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
    } else
        Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

    phase = 1 - phase;
    return Z;
}

void perturbVector(scs_float *v, scs_int l) {
    scs_int i;
    for (i = 0; i < l; i++) {
        v[i] += 0.01 * rand_gauss();
    }
}

void genRandomProbData(scs_int nnz, scs_int col_nnz, Data *d, Cone *k,
                       Sol *opt_sol) {
    scs_int n = d->n;
    scs_int m = d->m;
    AMatrix *A = d->A = scs_calloc(1, sizeof(AMatrix));
    scs_float *b = d->b = scs_calloc(m, sizeof(scs_float));
    scs_float *c = d->c = scs_calloc(n, sizeof(scs_float));
    scs_float *x = opt_sol->x = scs_calloc(n, sizeof(scs_float));
    scs_float *y = opt_sol->y = scs_calloc(m, sizeof(scs_float));
    scs_float *s = opt_sol->s = scs_calloc(m, sizeof(scs_float));
    /* temporary variables */
    scs_float *z = scs_calloc(m, sizeof(scs_float));
    scs_int i, j, r;

    A->i = scs_calloc(nnz, sizeof(scs_int));
    A->p = scs_calloc((n + 1), sizeof(scs_int));
    A->x = scs_calloc(nnz, sizeof(scs_float));
    A->n = d->n;
    A->m = d->m;
    /* y, s >= 0 and y'*s = 0 */
    for (i = 0; i < m; i++) {
        y[i] = z[i] = rand_scs_float();
    }

    projDualCone(y, k, SCS_NULL, SCS_NULL, -1);

    for (i = 0; i < m; i++) {
        b[i] = s[i] = y[i] - z[i];
    }

    for (i = 0; i < n; i++) {
        x[i] = rand_scs_float();
    }

    /* 	c = -A'*y
     b = A*x + s
     */
    A->p[0] = 0;
    scs_printf("Generating random matrix:\n");
    /*
     TODO: this only works probabilistically, ok for low density matrices
     */
    for (j = 0; j < n; j++) { /* column */
        if (j * 100 % n == 0 && (j * 100 / n) % 10 == 0) {
            scs_printf("%ld%%\n", (long)(j * 100 / n));
        }
        for (r = 0; r < col_nnz; r++) { /* row index */
            i = rand() % m;             /* row */
            A->x[r + j * col_nnz] = rand_scs_float();
            A->i[r + j * col_nnz] = i;

            b[i] += A->x[r + j * col_nnz] * x[j];

            c[j] -= A->x[r + j * col_nnz] * y[i];
        }
        A->p[j + 1] = (j + 1) * col_nnz;
    }
    scs_printf("done\n");
    scs_free(z);
}

#endif
