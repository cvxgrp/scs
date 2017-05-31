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
     */
    scs_float innerProd(const scs_float *x, const scs_float *y, scs_int len);

    /**
     * Returns the square Euclidean norm of a vector.
     * @param v
     * @param len
     * @return 
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
     * Returns the modified sign of \c x, that is
     * \f[
     *   \mathrm{sgn}(x) = 
     *  \cases{
     *     1, \text{ if }\, x \geq 0,\cr
     *     0, \text{ otherwise }\cr
     *  }
     * \f]
     * @param x
     * @return 
     */
    scs_float scs_sgn(scs_float x);

    /**
     * Returns the absolute value of \c x.
     * @param x
     * @return 
     */
    scs_float scs_abs(scs_float x);

#ifdef __cplusplus
}
#endif
#endif
