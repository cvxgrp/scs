#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

/* this struct defines the data matrix A */
struct A_DATA_MATRIX {
    /* A is supplied in column compressed format */
    scs_float *x; /* A values, size: NNZ A */
    scs_int *i;   /* A row index, size: NNZ A */
    scs_int *p;   /* A column pointer, size: n+1 */
    scs_int m, n; /* m rows, n cols */
};

#ifdef __cplusplus
}
#endif
#endif
