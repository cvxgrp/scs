#include "linAlg.h"
#include <math.h>

/* x = b*a */
void setAsScaledArray(scs_float *x, const scs_float *a, const scs_float b,
                      scs_int len) {
    scs_int i;
    for (i = 0; i < len; ++i)
        x[i] = b * a[i];
}

/* a *= b */
void scaleArray(scs_float *a, const scs_float b, scs_int len) {
    scs_int i;
    for (i = 0; i < len; ++i)
        a[i] *= b;
}

/* x'*y */
scs_float innerProd(const scs_float *x, const scs_float *y, scs_int len) {
    scs_int i;
    scs_float ip = 0.0;
    for (i = 0; i < len; ++i) {
        ip += x[i] * y[i];
    }
    return ip;
}

/* ||v||_2^2 */
scs_float calcNormSq(const scs_float *v, scs_int len) {
    scs_int i;
    scs_float nmsq = 0.0;
    for (i = 0; i < len; ++i) {
        nmsq += v[i] * v[i];
    }
    return nmsq;
}

/* ||v||_2 */
scs_float calcNorm(const scs_float *v, scs_int len) {
    return SQRTF(calcNormSq(v, len));
}

scs_float calcNormInf(const scs_float *a, scs_int l) {
    scs_float tmp, max = 0.0;
    scs_int i;
    for (i = 0; i < l; ++i) {
        tmp = ABS(a[i]);
        if (tmp > max)
            max = tmp;
    }
    return max;
}

/* saxpy a += sc*b */
void addScaledArray(scs_float *a, const scs_float *b, scs_int n,
                    const scs_float sc) {
    scs_int i;
    for (i = 0; i < n; ++i) {
        a[i] += sc * b[i];
    }
}

scs_float calcNormDiff(const scs_float *a, const scs_float *b, scs_int l) {
    scs_float nmDiff = 0.0, tmp;
    scs_int i;
    for (i = 0; i < l; ++i) {
        tmp = (a[i] - b[i]);
        nmDiff += tmp * tmp;
    }
    return SQRTF(nmDiff);
}

scs_float calcNormInfDiff(const scs_float *a, const scs_float *b, scs_int l) {
    scs_float tmp, max = 0.0;
    scs_int i;
    for (i = 0; i < l; ++i) {
        tmp = ABS(a[i] - b[i]);
        if (tmp > max)
            max = tmp;
    }
    return max;
}
