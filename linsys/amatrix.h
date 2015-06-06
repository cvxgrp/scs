#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

/* this struct defines the data matrix A */
struct A_DATA_MATRIX {
	/* A is supplied in dense format, COLUMN major order */
	scs_float * x;  /* A values, size: m by n */
    scs_int m, n;   /* m rows, n cols */
};

#ifdef __cplusplus
}
#endif
#endif
