#ifndef CS_H_GUARD
#define CS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

typedef struct SCS(cs_sparse) /* matrix in compressed-column or triplet form */
{
  scs_int nzmax; /* maximum number of entries */
  scs_int m;     /* number of rows */
  scs_int n;     /* number of columns */
  scs_int *p;    /* column pointers (size n+1) or col indices (size nzmax) */
  scs_int *i;    /* row indices, size nzmax */
  scs_float *x;  /* numerical values, size nzmax */
  scs_int nz;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs;

cs *SCS(cs_compress)(const cs *T);
cs *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                    scs_int triplet);
cs *SCS(cs_spfree)(cs *A);
scs_float SCS(cs_cumsum)(scs_int *p, scs_int *c, scs_int n);
scs_int *SCS(cs_pinv)(scs_int const *p, scs_int n);
cs *SCS(cs_symperm)(const cs *A, const scs_int *pinv, scs_int values);

#ifdef __cplusplus
}
#endif
#endif
