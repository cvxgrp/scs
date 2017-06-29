#ifndef CS_H_GUARD
#define CS_H_GUARD

#include "glbopts.h"

/**
 * \brief Matrix in compressed-column or triplet form.
 *
 */
typedef struct cs_sparse 
    {
    scs_int nzmax; /**< \brief maximum number of entries */
    scs_int m;     /**< \brief number of rows */
    scs_int n;     /**< \brief number of columns */
    scs_int *p;    /**< \brief column pointers (size n+1) or col indices (size nzmax) */
    scs_int *i;    /**< \brief row indices, size nzmax */
    scs_float *x;  /**< \brief numerical values, size nzmax */
    scs_int nz;    /**< \brief Number of entries in triplet matrix, -1 for compressed-col */
} cs;

/**
 * \brief Compress a triplet matrix into a column-packed representation.
 */
cs *cs_compress(const cs *T);
/**
 * \brief Frees the memory of <code>x</code> and <code>w></code>.
 * 
 * If <code>ok</code> is nonzero, it returns <code>C</code>, otherwise
 * it frees <code>C</code> (it calls ::cs_spfree) and returns ::SCS_NULL.
 * 
 * @param C
 * @param w
 * @param x
 * @param ok
 * @return 
 */
cs *cs_done(cs *C, void *w, void *x, scs_int ok);
/**
 * \brief Allocates a sparse matrix of given dimensions.
 * 
 * @param m
 * @param n
 * @param nzmax
 * @param values
 * @param triplet
 * @return 
 */
cs *cs_spalloc(scs_int m, scs_int n, scs_int nzmax, scs_int values,
               scs_int triplet);
cs *cs_spfree(cs *A);
scs_float cs_cumsum(scs_int *p, scs_int *c, scs_int n);
scs_int *cs_pinv(scs_int const *p, scs_int n);
cs *cs_symperm(const cs *A, const scs_int *pinv, scs_int values);
#endif
