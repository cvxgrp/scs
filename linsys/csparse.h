/* Routines modified from CSparse, T. Davis et al */

#ifndef CS_H_GUARD
#define CS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"

/* matrix in compressed-column or triplet form */
typedef struct SPARSE_MATRIX {
  scs_int nzmax; /* maximum number of entries */
  scs_int m;     /* number of rows */
  scs_int n;     /* number of columns */
  scs_int *p;    /* column pointers (size n+1) or col indices (size nzmax) */
  scs_int *i;    /* row indices, size nzmax */
  scs_float *x;  /* numerical values, size nzmax */
  scs_int nz;    /* # of entries in triplet matrix, -1 for compressed-col */
} csc;

csc *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                     scs_int triplet);
csc *SCS(cs_done)(csc *C, void *w, void *x, scs_int ok);
csc *SCS(cs_compress)(const csc *T, scs_int *idx_mapping);
scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);
csc *SCS(cs_spfree)(csc *A);

#ifdef __cplusplus
}
#endif
#endif
