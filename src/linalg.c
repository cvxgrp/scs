#include "linalg.h"
#include "scs_blas.h"
#include <math.h>

/* these routines do not have BLAS implementations (that I can find at least) */

scs_float SCS(norm_diff)(const scs_float *a, const scs_float *b, scs_int len) {
  scs_float nm_diff = 0.0, tmp;
  scs_int i;
  for (i = 0; i < len; ++i) {
    tmp = (a[i] - b[i]);
    nm_diff += tmp * tmp;
  }
  return SQRTF(nm_diff);
}

scs_float SCS(norm_inf_diff)(const scs_float *a, const scs_float *b,
                             scs_int len) {
  scs_float tmp, max = 0.0;
  scs_int i;
  for (i = 0; i < len; ++i) {
    tmp = ABS(a[i] - b[i]);
    if (tmp > max) {
      max = tmp;
    }
  }
  return max;
}

#ifndef USE_LAPACK
/* Self-rolled basic linear algebra routines */

/* a *= b */
void SCS(scale_array)(scs_float *a, const scs_float b, scs_int len) {
  scs_int i;
  for (i = 0; i < len; ++i)
    a[i] *= b;
}

/* x'*y */
scs_float SCS(dot)(const scs_float *x, const scs_float *y, scs_int len) {
  scs_int i;
  scs_float ip = 0.0;
  for (i = 0; i < len; ++i) {
    ip += x[i] * y[i];
  }
  return ip;
}

/* ||v||_2^2 */
scs_float SCS(norm_sq)(const scs_float *v, scs_int len) {
  scs_int i;
  scs_float nmsq = 0.0;
  for (i = 0; i < len; ++i) {
    nmsq += v[i] * v[i];
  }
  return nmsq;
}

/* ||v||_2 */
scs_float SCS(norm_2)(const scs_float *v, scs_int len) {
  return SQRTF(SCS(norm_sq)(v, len));
}

scs_float SCS(norm_inf)(const scs_float *a, scs_int len) {
  scs_float tmp, max = 0.0;
  scs_int i;
  for (i = 0; i < len; ++i) {
    tmp = ABS(a[i]);
    if (tmp > max) {
      max = tmp;
    }
  }
  return max;
}

/* axpy a += sc*b */
void SCS(add_scaled_array)(scs_float *a, const scs_float *b, scs_int n,
                           const scs_float sc) {
  scs_int i;
  for (i = 0; i < n; ++i) {
    a[i] += sc * b[i];
  }
}

scs_float SCS(mean)(const scs_float *x, scs_int n) {
  scs_int i;
  scs_float mean = 0.;
  for (i = 0; i < n; ++i) {
    mean += x[i];
  }
  return mean / n;
}

#else
/* If we have BLAS / LAPACK we may as well use them */

#ifdef __cplusplus
extern "C" {
#endif

scs_float BLAS(nrm2)(blas_int *n, const scs_float *x, blas_int *incx);
scs_float BLAS(dot)(const blas_int *n, const scs_float *x, const blas_int *incx,
                    const scs_float *y, const blas_int *incy);
scs_float BLAS(lange)(const char *norm, const blas_int *m, const blas_int *n,
                      const scs_float *a, blas_int *lda, scs_float *work);
void BLAS(axpy)(blas_int *n, const scs_float *a, const scs_float *x,
                blas_int *incx, scs_float *y, blas_int *incy);
void BLAS(scal)(const blas_int *n, const scs_float *sa, scs_float *sx,
                const blas_int *incx);

#ifdef __cplusplus
}
#endif

/* a *= b */
void SCS(scale_array)(scs_float *a, const scs_float b, scs_int len) {
  blas_int bone = 1;
  blas_int blen = (blas_int)len;
  BLAS(scal)(&blen, &b, a, &bone);
}

/* x'*y */
scs_float SCS(dot)(const scs_float *x, const scs_float *y, scs_int len) {
  blas_int bone = 1;
  blas_int blen = (blas_int)len;
  return BLAS(dot)(&blen, x, &bone, y, &bone);
}

/* ||v||_2^2 */
scs_float SCS(norm_sq)(const scs_float *v, scs_int len) {
  scs_float nrm = SCS(norm_2)(v, len);
  return nrm * nrm;
}

/* ||v||_2 */
scs_float SCS(norm_2)(const scs_float *v, scs_int len) {
  blas_int bone = 1;
  blas_int blen = (blas_int)len;
  return BLAS(nrm2)(&blen, v, &bone);
}

scs_float SCS(norm_inf)(const scs_float *a, scs_int len) {
  blas_int bone = 1;
  blas_int blen = (blas_int)len;
  return BLAS(lange)("Max", &blen, &bone, a, &bone, SCS_NULL);
}

/* axpy a += sc*b */
void SCS(add_scaled_array)(scs_float *a, const scs_float *b, scs_int len,
                           const scs_float sc) {
  blas_int bone = 1;
  blas_int blen = (blas_int)len;
  BLAS(axpy)(&blen, &sc, b, &bone, a, &bone);
}

scs_float SCS(mean)(const scs_float *x, scs_int n) {
  blas_int bone = 1;
  blas_int bzero = 0;
  blas_int blen = (blas_int)n;
  scs_float y = 1.0;
  return BLAS(dot)(&blen, x, &bone, &y, &bzero) / n;
}

#endif
