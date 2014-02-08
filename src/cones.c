#include "cones.h"
#ifdef LAPACK_LIB_FOUND
#include <cblas.h>
#include <lapacke.h>
#endif

/* private data to help cone projection step */
static struct ConeData_t {
	/* workspace for eigenvector decompositions: */
	pfloat * Xs, *Z, *e;
} c;

void projExpCone(pfloat * v);

idxint initCone(Cone * k) {
	if (k->ssize && k->s) {
		/* eigenvector decomp workspace */
		idxint i, nMax = 0;
		for (i = 0; i < k->ssize; ++i) {
			if (k->s[i] > nMax)
				nMax = k->s[i];
		}
		c.Xs = scs_calloc(nMax * nMax, sizeof(pfloat));
		c.Z = scs_calloc(nMax * nMax, sizeof(pfloat));
		c.e = scs_calloc(nMax, sizeof(pfloat));
		if (!c.Xs || !c.Z || !c.e)
			return -1;
	} else {
		c.Xs = NULL;
		c.Z = NULL;
		c.e = NULL;
	}
	return 0;
}

idxint validateCones(Cone * k) {
	idxint i;
	if (k->f && k->f < 0) {
		scs_printf("free cone error\n");
		return -1;
	}
	if (k->l && k->l < 0) {
		scs_printf("lp cone error\n");
		return -1;
	}
	if (k->qsize && k->q) {
		for (i = 0; i < k->qsize; ++i) {
			if (k->q[i] <= 0) {
				scs_printf("soc cone error\n");
				return -1;
			}
		}
	}
	if (k->ssize && k->s) {
		for (i = 0; i < k->ssize; ++i) {
			if (k->s[i] <= 0) {
				scs_printf("sd cone error\n");
				return -1;
			}
		}
	}
	if (k->ed && k->ed < 0) {
		scs_printf("ep cone error\n");
		return -1;
	}
	if (k->ep && k->ep < 0) {
		scs_printf("ed cone error\n");
		return -1;
	}
	return 0;
}

idxint getFullConeDims(Cone * k) {
	idxint i, c = 0;
	if (k->f)
		c += k->f;
	if (k->l)
		c += k->l;
	if (k->qsize && k->q) {
		for (i = 0; i < k->qsize; ++i) {
			c += k->q[i];
		}
	}
	if (k->ssize && k->s) {
		for (i = 0; i < k->ssize; ++i) {
			c += k->s[i] * k->s[i];
		}
	}
	if (k->ed)
		c += 3 * k->ed;
	if (k->ep)
		c += 3 * k->ep;
	return c;
}

void finishCone() {
	if (c.Xs)
		scs_free(c.Xs);
	if (c.Z)
		scs_free(c.Z);
	if (c.e)
		scs_free(c.e);
}

char * getConeHeader(Cone * k) {
	char * tmp = scs_malloc(sizeof(char) * 256);
	idxint i, socVars, socBlks, sdVars, sdBlks, expPvars, expDvars, len;
	socVars = 0;
	socBlks = 0;
	if (k->qsize && k->q) {
		socBlks = k->qsize;
		for (i = 0; i < k->qsize; i++) {
			socVars += k->q[i];
		}
	}
	sdVars = 0;
	sdBlks = 0;
	if (k->ssize && k->s) {
		sdBlks = k->ssize;
		for (i = 0; i < k->ssize; i++) {
			sdVars += k->s[i] * k->s[i];
		}
	}
	expPvars = 0;
	if (k->ep) {
		expPvars = 3 * k->ep;
	}
	expDvars = 0;
	if (k->ed) {
		expDvars = 3 * k->ed;
	}
	len = sprintf(tmp, "cones:\tzero/free vars: %i\n\tlinear vars: %i\n\tsoc vars: %i, soc blks: %i\n\tsd vars: %i, sd blks: %i\n\texp vars: %i, dual exp vars: %i\n",
					(int) (k->f ? k->f : 0), (int) (k->l ? k->l : 0), (int) socVars, (int) socBlks, (int) sdVars,
					(int) sdBlks, (int) expPvars, (int) expDvars);
	return tmp;
}

/* in place projection (with branches) */
void projCone(pfloat *x, Cone * k, idxint iter) {
	idxint i;
	idxint count = (k->f ? k->f : 0);

	if (k->l) {
		/* project onto positive orthant */
		for (i = count; i < count + k->l; ++i) {
			if (x[i] < 0.0)
				x[i] = 0.0;
			/*x[i] = (x[i] < 0.0) ? 0.0 : x[i]; */
		}
		count += k->l;
	}

	if (k->qsize && k->q) {
		/* project onto SOC */
		for (i = 0; i < k->qsize; ++i) {
			if (k->q[i] == 1) {
				if (x[count] < 0.0)
					x[count] = 0.0;
			} else {
				pfloat v1 = x[count];
				pfloat s = calcNorm(&(x[count + 1]), k->q[i] - 1);
				pfloat alpha = (s + v1) / 2.0;

				if (s <= v1) { /* do nothing */
				} else if (s <= -v1) {
					memset(&(x[count]), 0, k->q[i] * sizeof(pfloat));
				} else {
					x[count] = alpha;
					scaleArray(&(x[count + 1]), alpha / s, k->q[i] - 1);
					/*cblas_dscal(k->q[i]-1, alpha/s, &(x[count+1]),1); */
				}
			}
			count += k->q[i];
		}
	}

	if (k->ssize && k->s) {
#ifdef LAPACK_LIB_FOUND
		/* project onto PSD cone */
		for (i=0; i < k->ssize; ++i) {
			projectsdc(&(x[count]), k->s[i]);
			count += (k->s[i])*(k->s[i]);
		}
#else
		if (k->ssize > 0) {
			scs_printf("WARNING: solving SDP, but no blas/lapack libraries were linked!\n");
			scs_printf("scs will return nonsense!\n");
			for (i = 0; i < k->ssize; ++i) {
				scaleArray(&(x[count]), NAN, k->s[i] * k->s[i]);
				count += (k->s[i]) * (k->s[i]);
			}
		}
#endif
	}

	if (k->ep) {
		pfloat r, s, t;
		idxint idx;
		/*
		 * exponential cone is not self dual, if s \in K
		 * then y \in K^* and so if K is the primal cone
		 * here we project onto K^*, via Moreau
		 */
		scaleArray(&(x[count]), -1, 3 * k->ep); /* x = -x; */
		#pragma omp parallel for private(r,s,t,idx)
		for (i = 0; i < k->ep; ++i) {
			idx = count + 3 * i;
			r = x[idx];
			s = x[idx + 1];
			t = x[idx + 2];

			projExpCone(&(x[idx]));

			x[idx] -= r;
			x[idx + 1] -= s;
			x[idx + 2] -= t;
		}
		count += 3 * k->ep;
	}

	if (k->ed) {
		/* exponential cone: */
		#pragma omp parallel for
		for (i = 0; i < k->ed; ++i) {
			projExpCone(&(x[count + 3 * i]));
		}
		count += 3 * k->ed;
	}
	/* project onto OTHER cones */
}

pfloat expNewtonOneD(pfloat rho, pfloat y_hat, pfloat z_hat) {
	pfloat t = MAX(-z_hat, 1e-6);
	pfloat f, fp;
	idxint i;
	for (i = 0; i < 100; ++i) {

		f = t * (t + z_hat) / rho / rho - y_hat / rho + log(t / rho) + 1;
		fp = (2 * t + z_hat) / rho / rho + 1 / t;

		t = t - f / fp;

		if (t <= -z_hat) {
			return 0;
		} else if (t <= 0) {
			return z_hat;
		} else if ( ABS(f) < 1e-9) {
			break;
		}
	}
	return t + z_hat;
}

void expSolveForXWithRho(pfloat * v, pfloat * x, pfloat rho) {
	x[2] = expNewtonOneD(rho, v[1], v[2]);
	x[1] = (x[2] - v[2]) * x[2] / rho;
	x[0] = v[0] - rho;
}

pfloat expCalcGrad(pfloat * v, pfloat * x, pfloat rho) {
	expSolveForXWithRho(v, x, rho);
	if (x[1] <= 1e-12) {
		return x[0];
	} else {
		return x[0] + x[1] * log(x[1] / x[2]);
	}
}

void expGetRhoUb(pfloat * v, pfloat * x, pfloat * ub, pfloat * lb) {
	*lb = 0;
	*ub = 0.125;
	while (expCalcGrad(v, x, *ub) > 0) {
		*lb = *ub;
		(*ub) *= 2;
	}
}

/* project onto the exponential cone, v has dimension *exactly* 3 */
void projExpCone(pfloat * v) {
	idxint i;
	pfloat ub, lb, rho, g, x[3];
	pfloat r = v[0], s = v[1], t = v[2];

	/* v in cl(Kexp) */
	if ((s * exp(r / s) <= t && s > 0) || (r <= 0 && s == 0 && t >= 0)) {
		return;
	}

	/* -v in Kexp^* */
	if ((-r < 0 && r * exp(s / r) <= -exp(1) * t) || (-r == 0 && -s >= 0 && -t >= 0)) {
		memset(v, 0, 3 * sizeof(pfloat));
		return;
	}

	/* special case with analytical solution */
	if (r < 0 && s < 0) {
		v[1] = 0.0;
		v[2] = MAX(v[2], 0);
		return;
	}
	expGetRhoUb(v, x, &ub, &lb);
	for (i = 0; i < 100; ++i) {
		rho = (ub + lb) / 2;
		g = expCalcGrad(v, x, rho);
		if (g > 0) {
			lb = rho;
		} else {
			ub = rho;
		}
		if (ub - lb < 1e-9) {
			break;
		}
	}
	v[0] = x[0];
	v[1] = x[1];
	v[2] = x[2];
}

#ifdef LAPACK_LIB_FOUND
void projectsdc(pfloat *X, idxint n)
{ /* project onto the positive semi-definite cone */
	if (n == 1) {
		if(X[0] < 0.0) X[0] = 0.0;
		return;
	}

	idxint i, j;
	int m;
	pfloat * Xs = c.Xs;
	pfloat * Z = c.Z;
	pfloat * e = c.e;
	memcpy(Xs,X,n*n*sizeof(pfloat));

	/* Xs = X + X', save div by 2 for eigen-recomp */
	for (i = 0; i < n; ++i) {
		cblas_daxpy(n, 1, &(X[i]), n, &(Xs[i*n]), 1);
		/*b_daxpy(n, 1, &(X[i]), n, &(Xs[i*n]), 1); */
	}

	pfloat EIG_TOL = 1e-8;
	pfloat vupper = calcNorm(Xs,n*n);
	m = 0;
	LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'V', 'U', n, Xs, n, 0.0, vupper, -1, -1, EIG_TOL, &m, e, Z, n , NULL);
	memset(X, 0, n*n*sizeof(pfloat));
	for (i = 0; i < m; ++i) {
		cblas_dsyr(CblasColMajor, CblasLower, n, e[i]/2, &(Z[i*n]), 1, X, n);
		/* b_dsyr('L', n, -e[i]/2, &(Z[i*n]), 1, Xs, n); */
	}
	/* fill in upper half  */
	for (i = 0; i < n; ++i) {
		for (j = i+1; j < n; ++j) {
			X[i + j*n] = X[j + i*n];
		}
	}
}
#endif

