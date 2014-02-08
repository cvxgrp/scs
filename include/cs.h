#ifndef CS_H_GUARD
#define CS_H_GUARD

#include "glbopts.h"

typedef struct cs_sparse /* matrix in compressed-column or triplet form */
{
	idxint nzmax; /* maximum number of entries */
	idxint m; /* number of rows */
	idxint n; /* number of columns */
	idxint *p; /* column pointers (size n+1) or col indices (size nzmax) */
	idxint *i; /* row indices, size nzmax */
	pfloat *x; /* numerical values, size nzmax */
	idxint nz; /* # of entries in triplet matrix, -1 for compressed-col */
} cs;

cs *cs_compress(const cs *T);
cs *cs_done(cs *C, void *w, void *x, idxint ok);
cs *cs_spalloc(idxint m, idxint n, idxint nzmax, idxint values, idxint triplet);
cs *cs_spfree(cs *A);
pfloat cs_cumsum(idxint *p, idxint *c, idxint n);
idxint *cs_pinv(idxint const *p, idxint n);
cs *cs_symperm(const cs *A, const idxint *pinv, idxint values);
#endif
