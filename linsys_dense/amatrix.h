#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

/* this struct defines the data matrix A */
struct A_DATA_MATRIX {
	/* A is supplied in dense format, COLUMN major order */
	pfloat * x; /* A values, size: n*m */
};

#endif
