#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include <math.h>

    /**
     * Performs the operation
     * \f[
     *  x \leftarrow b\cdot a,
     * \f]
     * where \c a is a vector and \c b is a scalar.
     * 
     * @param x
     * @param a
     * @param b
     * @param len
     * 
     * \note with loop unrolling for speed
     */
    void setAsScaledArray(
            scs_float *x,
            const scs_float *a,
            const scs_float b,
            scs_int len);

    /**
     * Performs the operation
     * \f[
     *   a \leftarrow b\cdot a
     * \f]
     * @param a
     * @param b
     * @param len
     * 
     * \note with loop unrolling for speed
     */
    void scaleArray(scs_float *a, const scs_float b, scs_int len);

    /**
     * Computes the inner product of two vectors, that is
     * \f[
     *  \langle x, y \rangle = x'y = \sum_{i=1}^{\mathrm{len}}x_i y_i.
     * \f]
     * @param x
     * @param y
     * @param len
     * @return 
     * 
     * \note with loop unrolling for speed
     */
    scs_float innerProd(const scs_float *x, const scs_float *y, scs_int len);

    /**
     * Returns the square Euclidean norm of a vector.
     * @param v
     * @param len
     * @return 
     * 
     * \note uses ::innerProd
     */
    scs_float calcNormSq(const scs_float *v, scs_int len);

    /**
     * Returns the Euclidean norm of a vector.
     * @param v
     * @param len
     * @return 
     */
    scs_float calcNorm(const scs_float *v, scs_int len);

    /**
     * Returns the infinity norm of a vector.
     * @param a
     * @param l
     * @return 
     */
    scs_float calcNormInf(const scs_float *a, scs_int l);


    /**
     * Performs the operation
     * \f[
     *  a \leftarrow a + \gamma b
     * \f]
     * @param a vector \c a
     * @param b vector \c b
     * @param n length of \c a
     * @param sc the scalar \f$\gamma\f$
     * 
     * \note with loop unrolling for speed
     */
    void addScaledArray(
            scs_float *a,
            const scs_float *b,
            scs_int n,
            const scs_float sc);

    /**
     * Performs the operation
     * \f[
     *  a \leftarrow a + b
     * \f]
     * @param a vector \c a
     * @param b vector \c b
     * @param n length of \c a
     * 
     * \note with loop unrolling for speed
     */
    void addArray(scs_float *a, const scs_float *b, scs_int n);

    /**
     * Performs the operation
     * \f[
     *  a \leftarrow a - b
     * \f]
     * @param a vector \c a
     * @param b vector \c b
     * @param n length of \c a
     * 
     * \note with loop unrolling for speed
     */
    void subtractArray(scs_float *a, const scs_float *b, scs_int n);

    /**
     * Returns the Euclidean norm of the difference of two vectors
     * @param a
     * @param b
     * @param l
     * @return 
     */
    scs_float calcNormDiff(const scs_float *a, const scs_float *b, scs_int l);

    /**
     * Returns the infinity norm of the difference of two vectors
     * @param a
     * @param b
     * @param l
     * @return 
     */
    scs_float calcNormInfDiff(const scs_float *a, const scs_float *b, scs_int l);
    

    /**
     * Perofrms the operation
     * \f[
     *    C \leftarrow \beta C + \alpha A B
     * \f]
     * 
     * @param m number of rows of matrix \f$A\f$
     * @param n number of columns of matrix \f$B\f$
     * @param k number of rows of matrix \f$B\f$ (columns of \f$A\f$)
     * @param alpha coefficient \f$\alpha\f$
     * @param A pointer to matrix \f$A\f$
     * @param incRowA increment in traversing the rows of \f$A\f$
     * @param incColA increment in traversing the columns of \f$A\f$
     * @param B pointer to matrix \f$B\f$
     * @param incRowB increment in traversing the rows of \f$B\f$
     * @param incColB increment in traversing the columns of \f$B\f$
     * @param beta coefficient \f$\beta\f$
     * @param C pointer to matrix \f$C\f$
     * @param incRowC increment in traversing the rows of \f$C\f$
     * @param incColC increment in traversing the columns of \f$C\f$
     * 
     * @see ::matrixMultiplicationColumnPacked
     * 
     * \note The implementation of this method is that of 
     * [ULMBLAS](http://apfel.mathematik.uni-ulm.de/~lehn/sghpc/gemm/page13/index.html).
     * 
     * \note The original source code is available at 
     * [this link](http://apfel.mathematik.uni-ulm.de/~lehn/sghpc/gemm/page13/index.html).
     * 
     * \note The [ULMBLAS project](https://github.com/michael-lehn/ulmBLAS) is available
     * on github and is licensed with the 
     * [new BSD licence](https://github.com/michael-lehn/ulmBLAS/blob/master/LICENSE).
     * 
     * \warning This function works only with \c double precision data.
     * 
     */
    void dgemm_nn(
            int m,
            int n,
            int k,
            double alpha,
            const double *A,
            int incRowA,
            int incColA,
            const double *B,
            int incRowB,
            int incColB,
            double beta,
            double *C,
            int incRowC,
            int incColC);

    /**
     * Perofrms the operation \f$C \leftarrow \beta C + \alpha A B,\f$
     * where \f$A\f$, \f$B\f$ and \f$C\f$ are column-packed matrices.
     * 
     * This method is a proxy to ::dgemm_nn.
     * 
     * @param m number of rows of matrix \f$A\f$
     * @param n number of columns of matrix \f$B\f$
     * @param k number of rows of matrix \f$B\f$ (columns of \f$A\f$)
     * @param alpha coefficient \f$\alpha\f$
     * @param A pointer to matrix \f$A\f$ in column-packed form
     * @param beta coefficient \f$\beta\f$
     * @param B pointer to matrix \f$B\f$ in column-packed form
     * @param C pointer to matrix \f$C\f$ in column-packed form
     * 
     * @see ::dgemm_nn
     */
    void matrixMultiplicationColumnPacked(
            int m,
            int n,
            int k,
            double alpha,
            const double *A,
            double beta,
            const double *B,
            double *C);
    
    /**
     * Perofrms the operation \f$C \leftarrow \beta C + \alpha A^{\top} B,\f$
     * where \f$A\f$, \f$B\f$ and \f$C\f$ are column-packed matrices.
     * 
     * This method is a proxy to ::dgemm_nn.
     * 
     * @param m number of rows of matrix \f$A\f$
     * @param n number of columns of matrix \f$B\f$
     * @param k number of rows of matrix \f$B\f$ (columns of \f$A\f$)
     * @param alpha coefficient \f$\alpha\f$
     * @param A pointer to matrix \f$A\f$ in column-packed form
     * @param beta coefficient \f$\beta\f$
     * @param B pointer to matrix \f$B\f$ in column-packed form
     * @param C pointer to matrix \f$C\f$ in column-packed form
     * 
     * @see ::dgemm_nn
     */
    void matrixMultiplicationTransColumnPacked(
            int m,
            int n,
            int k,
            double alpha,
            const double *A,
            double beta,
            const double *B,
            double *C) ;

#ifdef __cplusplus
}
#endif
#endif
