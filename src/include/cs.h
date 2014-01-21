#ifndef CS_H_GUARD
#define CS_H_GUARD

//#include <string.h>
//#include <stdlib.h>

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_compress (const cs *T);
cs *cs_done (cs *C, void *w, void *x, int ok);
cs *cs_spalloc (int m, int n, int nzmax, int values, int triplet);
cs *cs_spfree (cs *A);
double cs_cumsum (int *p, int *c, int n);
int *cs_pinv (int const *p, int n);
cs *cs_symperm (const cs *A, const int *pinv, int values);
#endif
