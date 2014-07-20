#include "linAlg.h"
#include <math.h>

/* x = b*a */
void setAsScaledArray(pfloat *x, const pfloat * a, const pfloat b, idxint len) {
	idxint i;
	for (i = 0; i < len; ++i)
		x[i] = b * a[i];
}

/* a *= b */
void scaleArray(pfloat * a, const pfloat b, idxint len) {
	idxint i;
	for (i = 0; i < len; ++i)
		a[i] *= b;
}

/* x'*y */
pfloat innerProd(const pfloat * x, const pfloat * y, idxint len) {
	idxint i;
	pfloat ip = 0.0;
	for (i = 0; i < len; ++i) {
		ip += x[i] * y[i];
	}
	return ip;
}

/* ||v||_2^2 */
pfloat calcNormSq(const pfloat * v, idxint len) {
	idxint i;
	pfloat nmsq = 0.0;
	for (i = 0; i < len; ++i) {
		nmsq += v[i] * v[i];
	}
	return nmsq;
}

/* ||v||_2 */
pfloat calcNorm(const pfloat * v, idxint len) {
	return SQRTF(calcNormSq(v, len));
}

pfloat calcNormInf(const pfloat *a, idxint l) {
	pfloat tmp, max = 0.0;
	idxint i;
	for (i = 0; i < l; ++i) {
		tmp = ABS(a[i]);
		if (tmp > max)
			max = tmp;
	}
	return max;
}

/* saxpy a += sc*b */
void addScaledArray(pfloat * a, const pfloat * b, idxint n, const pfloat sc) {
	idxint i;
	for (i = 0; i < n; ++i) {
		a[i] += sc * b[i];
	}
}

pfloat calcNormDiff(const pfloat *a, const pfloat *b, idxint l) {
	pfloat nmDiff = 0.0, tmp;
	idxint i;
	for (i = 0; i < l; ++i) {
		tmp = (a[i] - b[i]);
		nmDiff += tmp * tmp;
	}
	return SQRTF(nmDiff);
}

pfloat calcNormInfDiff(const pfloat *a, const pfloat *b, idxint l) {
	pfloat tmp, max = 0.0;
	idxint i;
	for (i = 0; i < l; ++i) {
		tmp = ABS(a[i] - b[i]);
		if (tmp > max)
			max = tmp;
	}
	return max;
}
