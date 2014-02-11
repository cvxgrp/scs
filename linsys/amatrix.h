#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

/* this struct defines the data matrix A */
struct A_DATA_MATRIX {
	/* A is supplied in column compressed format */
	pfloat * x; /* A values, size: NNZ A */
	idxint * i; /* A row index, size: NNZ A */
	idxint * p; /* A column pointer, size: n+1 */
};

#endif
