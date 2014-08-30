#include "cones.h"

#define CONE_RATE (2)
#define CONE_TOL (1e-7)
#define EXP_CONE_MAX_ITERS (100)

#ifdef LAPACK_LIB_FOUND
/* underscore for blas / lapack, single or double precision */
#ifdef NOBLASUNDERSCORE
#ifndef FLOAT
#define BLAS(x) d ## x
#else
#define BLAS(x) s ## x
#endif
#else
#ifndef FLOAT
#define BLAS(x) d ## x ## _
#else
#define BLAS(x) s ## x ## _
#endif
#endif

#ifdef MATLAB_MEX_FILE
typedef ptrdiff_t blasint;
#elif defined BLAS64
#include <stdint.h>
typedef int64_t blasint;
#else
typedef int blasint;
#endif

void BLAS(syevr)(char* jobz, char* range, char* uplo, blasint* n, scs_float* a, blasint* lda, scs_float* vl,
		scs_float* vu, blasint* il, blasint* iu, scs_float* abstol, blasint* m, scs_float* w, scs_float* z, blasint* ldz,
		blasint* isuppz, scs_float* work, blasint* lwork, blasint* iwork, blasint* liwork, blasint* info);
void BLAS(syr)(const char *uplo, const blasint *n, const scs_float *alpha, const scs_float *x, const blasint *incx,
		scs_float *a, const blasint *lda);
void BLAS(axpy)(const blasint *n, const scs_float *alpha, const scs_float *dx, const blasint *incx, scs_float *dy,
		const blasint *incy);
/* private data to help cone projection step */
static struct ConeData_t {
	/* workspace for eigenvector decompositions: */
	scs_float * Xs, *Z, *e, *work;
	blasint *iwork, lwork, liwork;
}c;
#endif

static timer coneTimer;
static scs_float totalConeTime;

 /*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size strictly larger than 1
 * returns length of boundaries array, boundaries malloc-ed here so should be freed
 */
scs_int getConeBoundaries(Cone * k, scs_int ** boundaries) {
	scs_int i, count = 0;
	scs_int len = 1 + k->qsize + k->ssize + k->ed + k->ep;
	scs_int * b = scs_malloc(sizeof(scs_int) * len);
	b[count] = k->f + k->l;
	count += 1;
	if (k->qsize > 0) {
		memcpy(&b[count], k->q, k->qsize * sizeof(scs_int));
	}
	count += k->qsize;
	for (i = 0; i < k->ssize; ++i) {
		b[count + i] = k->s[i] * k->s[i];
	}
	count += k->ssize;
	for (i = 0; i < k->ep + k->ed; ++i) {
		b[count + i] = 3;
	}
	count += k->ep + k->ed;
	*boundaries = b;
	return len;
}

scs_int getFullConeDims(Cone * k) {
	scs_int i, c = 0;
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

scs_int validateCones(Data * d, Cone * k) {
	scs_int i;
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
			if (k->q[i] < 0) {
				scs_printf("soc cone error\n");
				return -1;
			}
		}
	}
	if (k->ssize && k->s) {
		for (i = 0; i < k->ssize; ++i) {
			if (k->s[i] < 0) {
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
	if (getFullConeDims(k) != d->m) {
		scs_printf("cone dimensions %i not equal to num rows in A = m = %i\n", (int) getFullConeDims(k), (int) d->m);
		return -1;
	}
	return 0;
}

char * getConeSummary(Info * info) {
	char * str = scs_malloc(sizeof(char) * 64);
	sprintf(str, "\tCones: avg projection time: %1.2es\n", totalConeTime / (info->iter + 1) / 1e3);
	totalConeTime = 0.0;
	return str;
}

void finishCone() {
#ifdef LAPACK_LIB_FOUND
	if (c.Xs)
	scs_free(c.Xs);
	if (c.Z)
	scs_free(c.Z);
	if (c.e)
	scs_free(c.e);
	if (c.work)
	scs_free(c.work);
	if (c.iwork)
	scs_free(c.iwork);
#endif
}

char * getConeHeader(Cone * k) {
	char * tmp = scs_malloc(sizeof(char) * 512);
	scs_int i, socVars, socBlks, sdVars, sdBlks, expPvars, expDvars;
    sprintf(tmp, "Cones:");
    if (k->f) {
        sprintf(tmp + strlen(tmp), "\tprimal zero / dual free vars: %i\n", (int) k->f);
    }
    if (k->l) {
        sprintf(tmp + strlen(tmp), "\tlinear vars: %i\n", (int) k->l);
    }
	socVars = 0;
	socBlks = 0;
	if (k->qsize && k->q) {
		socBlks = k->qsize;
		for (i = 0; i < k->qsize; i++) {
			socVars += k->q[i];
		}
        sprintf(tmp + strlen(tmp), "\tsoc vars: %i, soc blks: %i\n", (int) socVars, (int) socBlks);
	}
	sdVars = 0;
	sdBlks = 0;
	if (k->ssize && k->s) {
		sdBlks = k->ssize;
		for (i = 0; i < k->ssize; i++) {
			sdVars += k->s[i] * k->s[i];
		}
        sprintf(tmp + strlen(tmp), "\tsd vars: %i, sd blks: %i\n", (int) sdVars, (int) sdBlks);
	}
    if (k->ep || k->ed) {
	    expPvars = k->ep ? 3 * k->ep : 0;
	    expDvars = k->ed ? 3 * k->ed : 0;
        sprintf(tmp + strlen(tmp), "\texp vars: %i, dual exp vars: %i\n", (int) expPvars, (int) expDvars);
    }
	return tmp;
}

scs_int isSimpleSemiDefiniteCone(scs_int * s, scs_int ssize) {
	scs_int i;
	for (i = 0; i < ssize; i++) {
		if (s[i] >= 3) {
			return 0; /* false */
		}
	}
	return 1; /* true */
}

scs_float expNewtonOneD(scs_float rho, scs_float y_hat, scs_float z_hat) {
	scs_float t = MAX(-z_hat, 1e-6);
	scs_float f, fp;
	scs_int i;
	for (i = 0; i < EXP_CONE_MAX_ITERS; ++i) {

		f = t * (t + z_hat) / rho / rho - y_hat / rho + log(t / rho) + 1;
		fp = (2 * t + z_hat) / rho / rho + 1 / t;

		t = t - f / fp;

		if (t <= -z_hat) {
			return 0;
		} else if (t <= 0) {
			return z_hat;
		} else if ( ABS(f) < CONE_TOL) {
			break;
		}
	}
	return t + z_hat;
}

void expSolveForXWithRho(scs_float * v, scs_float * x, scs_float rho) {
	x[2] = expNewtonOneD(rho, v[1], v[2]);
	x[1] = (x[2] - v[2]) * x[2] / rho;
	x[0] = v[0] - rho;
}

scs_float expCalcGrad(scs_float * v, scs_float * x, scs_float rho) {
	expSolveForXWithRho(v, x, rho);
	if (x[1] <= 1e-12) {
		return x[0];
	} else {
		return x[0] + x[1] * log(x[1] / x[2]);
	}
}

void expGetRhoUb(scs_float * v, scs_float * x, scs_float * ub, scs_float * lb) {
	*lb = 0;
	*ub = 0.125;
	while (expCalcGrad(v, x, *ub) > 0) {
		*lb = *ub;
		(*ub) *= 2;
	}
}

/* project onto the exponential cone, v has dimension *exactly* 3 */
static scs_int projExpCone(scs_float * v, scs_int iter) {
	scs_int i;
	scs_float ub, lb, rho, g, x[3];
	scs_float r = v[0], s = v[1], t = v[2];
	scs_float tol = CONE_TOL; /* iter < 0 ? CONE_TOL : MAX(CONE_TOL, 1 / POWF((iter + 1), CONE_RATE)); */

	/* v in cl(Kexp) */
	if ((s * exp(r / s) <= t && s > 0) || (r <= 0 && s == 0 && t >= 0)) {
		return 0;
	}

	/* -v in Kexp^* */
	if ((-r < 0 && r * exp(s / r) <= -exp(1) * t) || (-r == 0 && -s >= 0 && -t >= 0)) {
		memset(v, 0, 3 * sizeof(scs_float));
		return 0;
	}

	/* special case with analytical solution */
	if (r < 0 && s < 0) {
		v[1] = 0.0;
		v[2] = MAX(v[2], 0);
		return 0;
	}

    /* iterative procedure to find projection, bisects on dual variable: */
	expGetRhoUb(v, x, &ub, &lb); /* get starting upper and lower bounds */
	for (i = 0; i < EXP_CONE_MAX_ITERS; ++i) {
		rho = (ub + lb) / 2; /* halfway between upper and lower bounds */
		g = expCalcGrad(v, x, rho); /* calculates gradient wrt dual var */
		if (g > 0) {
			lb = rho;
		} else {
			ub = rho;
		}
		if (ub - lb < tol) {
			break;
		}
	}
	/*
#ifdef EXTRAVERBOSE
	scs_printf("exponential cone proj iters %i\n", i);
#endif
	 */
	v[0] = x[0];
	v[1] = x[1];
	v[2] = x[2];
	return 0;
}

scs_int initCone(Cone * k) {
#ifdef LAPACK_LIB_FOUND
	scs_int i;
	blasint nMax = 0;
	scs_float eigTol = 1e-8;
	blasint negOne = -1;
	blasint m = 0;
    blasint info;
	scs_float wkopt;
	c.Xs = NULL;
	c.Z = NULL;
	c.e = NULL;
	c.work = NULL;
	c.iwork = NULL;
#endif
	totalConeTime = 0.0;
#ifdef EXTRAVERBOSE
    scs_printf("initCone\n");
#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
#endif

if (k->ssize && k->s) {
		if (isSimpleSemiDefiniteCone(k->s, k->ssize)) {
			return 0;
		}
#ifdef LAPACK_LIB_FOUND
		/* eigenvector decomp workspace */
		for (i = 0; i < k->ssize; ++i) {
			if (k->s[i] > nMax) {
				nMax = (blasint) k->s[i];
			}
		}
		c.Xs = scs_calloc(nMax * nMax, sizeof(scs_float));
		c.Z = scs_calloc(nMax * nMax, sizeof(scs_float));
		c.e = scs_calloc(nMax, sizeof(scs_float));

        BLAS(syevr)("Vectors", "All", "Upper", &nMax, c.Xs, &nMax, NULL, NULL, NULL, NULL,
            &eigTol, &m, c.e, c.Z, &nMax, NULL, &wkopt, &negOne, &(c.liwork), &negOne, &info);

        if (info != 0) {
            scs_printf("FATAL: syevr failure, info = %i\n", (int) info);
            return -1;
        }
        c.lwork = (blasint) (wkopt + 0.01); /* 0.01 for int casting safety */
        c.work = scs_malloc(c.lwork * sizeof(scs_float));
        c.iwork = scs_malloc(c.liwork * sizeof(blasint));

		if (!c.Xs || !c.Z || !c.e || !c.work || !c.iwork) {
			return -1;
		}
#else
		scs_printf("FATAL: Cannot solve SDPs with > 2x2 matrices without linked blas+lapack libraries\n");
		scs_printf("Edit scs.mk to point to blas+lapack libray locations\n");
		return -1;
#endif
	}
#ifdef EXTRAVERBOSE
    scs_printf("initCone complete\n");
#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
#endif
	return 0;
}

scs_int project2By2Sdc(scs_float *X) {
	scs_float a, b, d, l1, l2, x1, x2, rad;
	a = X[0];
	b = 0.5 * (X[1] + X[2]);
	d = X[3];

	rad = SQRTF((a - d) * (a - d) + 4 * b * b);
	/* l1 >= l2 always, since rad >= 0 */
	l1 = 0.5 * (a + d + rad);
	l2 = 0.5 * (a + d - rad);

	if (l2 >= 0) { /* both positive, just symmetrize */
		X[1] = b;
		X[2] = b;
		return 0;
	}
	if (l1 <= 0) { /* both negative, set to 0 */
		X[0] = 0;
		X[1] = 0;
		X[2] = 0;
		X[3] = 0;
		return 0;
	}
	/* l1 pos, l2 neg */
	x1 = 1 / SQRTF(1 + (l1 - a) * (l1 - a) / b / b);
	x2 = x1 * (l1 - a) / b;

	X[0] = l1 * x1 * x1;
	X[1] = l1 * x1 * x2;
	X[2] = X[1];
	X[3] = l1 * x2 * x2;
	return 0;
}

static scs_int projSemiDefiniteCone(scs_float *X, scs_int n, scs_int iter) {
	/* project onto the positive semi-definite cone */
#ifdef LAPACK_LIB_FOUND
	scs_int i, j;
	blasint one = 1;
	blasint m = 0;
	blasint nb = (blasint) n;
	scs_float * Xs = c.Xs;
	scs_float * Z = c.Z;
	scs_float * e = c.e;
	scs_float * work = c.work;
	blasint * iwork = c.iwork;
	blasint lwork = c.lwork;
	blasint liwork = c.liwork;

	scs_float eigTol = CONE_TOL; /* iter < 0 ? CONE_TOL : MAX(CONE_TOL, 1 / POWF(iter + 1, CONE_RATE)); */
	scs_float onef = 1.0;
	scs_float zero = 0.0;
	blasint info;
	scs_float vupper;
#endif
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		if (X[0] < 0.0) {
			X[0] = 0.0;
		}
		return 0;
	}
	if (n == 2) {
		return project2By2Sdc(X);
	}
#ifdef LAPACK_LIB_FOUND
	memcpy(Xs, X, n * n * sizeof(scs_float));

    /* Xs = X + X', save div by 2 for eigen-recomp */
    for (i = 0; i < n; ++i) {
        BLAS(axpy)(&nb, &onef, &(X[i]), &nb, &(Xs[i * n]), &one);
    }
    vupper = MAX(calcNorm(Xs, n * n), 0.001);
    /* Solve eigenproblem, reuse workspaces */
    BLAS(syevr)("Vectors", "VInterval", "Upper", &nb, Xs, &nb, &zero, &vupper,
            NULL, NULL, &eigTol, &m, e, Z, &nb, NULL, work, &lwork, iwork, &liwork, &info);
    if (info != 0) {
#ifdef EXTRAVERBOSE
        scs_printf("WARN: LAPACK syevr error, info = %i\n", info);
        scs_printf("syevr input parameter dump:\n");
        scs_printf("nb = %li\n", (long) nb);
        scs_printf("lwork = %li\n", (long) lwork);
        scs_printf("liwork = %li\n", (long) liwork);
        scs_printf("vupper = %f\n", vupper);
        scs_printf("eigTol = %e\n", eigTol);
        printArray(Xs, n * n, "Xs");
        printArray(X, n * n, "X");
        printArray(e, m, "e");
        printArray(Z, m * n, "Z");
#endif
    if (info < 0) return -1;
    }

	memset(X, 0, n * n * sizeof(scs_float));
	for (i = 0; i < m; ++i) {
		scs_float a = e[i] / 2;
		BLAS(syr)("Lower", &nb, &a, &(Z[i * n]), &one, X, &nb);
	}
	/* fill in upper half */
	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			X[i + j * n] = X[j + i * n];
		}
	}
#else
	scs_printf("FAILURE: solving SDP with > 2x2 matrices, but no blas/lapack libraries were linked!\n");
	scs_printf("scs will return nonsense!\n");
	scaleArray(X, NAN, n);
	return -1;
#endif
	return 0;
}

/* outward facing cone projection routine, iter is outer algorithm iteration, if iter < 0 then iter is ignored
    warm_start contains guess of projection (can be set to NULL) */
scs_int projDualCone(scs_float *x, Cone * k, const scs_float * warm_start, scs_int iter)  {
	scs_int i;
	scs_int count = (k->f ? k->f : 0);
#ifdef EXTRAVERBOSE
	timer projTimer;
	tic(&projTimer);
#endif
	tic(&coneTimer);


	if (k->l) {
		/* project onto positive orthant */
		for (i = count; i < count + k->l; ++i) {
			if (x[i] < 0.0)
				x[i] = 0.0;
			/*x[i] = (x[i] < 0.0) ? 0.0 : x[i]; */
		}
		count += k->l;
#ifdef EXTRAVERBOSE
	scs_printf("pos orthant proj time: %1.2es\n", tocq(&projTimer) / 1e3);
	tic(&projTimer);
#endif
	}

	if (k->qsize && k->q) {
		/* project onto SOC */
		for (i = 0; i < k->qsize; ++i) {
			if (k->q[i] == 0) {
				continue;
			}
			if (k->q[i] == 1) {
				if (x[count] < 0.0)
					x[count] = 0.0;
			} else {
				scs_float v1 = x[count];
				scs_float s = calcNorm(&(x[count + 1]), k->q[i] - 1);
				scs_float alpha = (s + v1) / 2.0;

				if (s <= v1) { /* do nothing */
				} else if (s <= -v1) {
					memset(&(x[count]), 0, k->q[i] * sizeof(scs_float));
				} else {
					x[count] = alpha;
					scaleArray(&(x[count + 1]), alpha / s, k->q[i] - 1);
				}
			}
			count += k->q[i];
		}
#ifdef EXTRAVERBOSE
	scs_printf("SOC proj time: %1.2es\n", tocq(&projTimer) / 1e3);
	tic(&projTimer);
#endif
	}

	if (k->ssize && k->s) {
		/* project onto PSD cone */
		for (i = 0; i < k->ssize; ++i) {
			if (k->s[i] == 0) {
				continue;
			}
			if (projSemiDefiniteCone(&(x[count]), k->s[i], iter) < 0) return -1;
			count += (k->s[i]) * (k->s[i]);
		}
#ifdef EXTRAVERBOSE
	scs_printf("SD proj time: %1.2es\n", tocq(&projTimer) / 1e3);
	tic(&projTimer);
#endif
	}

	if (k->ep) {
		scs_float r, s, t;
		scs_int idx;
		/*
		 * exponential cone is not self dual, if s \in K
		 * then y \in K^* and so if K is the primal cone
		 * here we project onto K^*, via Moreau
		 * \Pi_C^*(y) = y + \Pi_C(-y)
		 */
		scaleArray(&(x[count]), -1, 3 * k->ep); /* x = -x; */
#ifdef OPENMP
#pragma omp parallel for private(r,s,t,idx)
#endif
		for (i = 0; i < k->ep; ++i) {
			idx = count + 3 * i;
			r = x[idx];
			s = x[idx + 1];
			t = x[idx + 2];

			if (projExpCone(&(x[idx]), iter) < 0) return -1;

			x[idx] -= r;
			x[idx + 1] -= s;
			x[idx + 2] -= t;
		}
		count += 3 * k->ep;
#ifdef EXTRAVERBOSE
	scs_printf("EP proj time: %1.2es\n", tocq(&projTimer) / 1e3);
	tic(&projTimer);
#endif
	}


	if (k->ed) {
		/* exponential cone: */
#ifdef OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < k->ed; ++i) {
			if (projExpCone(&(x[count + 3 * i]), iter) < 0) return -1;
		}
		count += 3 * k->ed;
#ifdef EXTRAVERBOSE
	scs_printf("ED proj time: %1.2es\n", tocq(&projTimer) / 1e3);
	tic(&projTimer);
#endif
	}
	/* project onto OTHER cones */
	totalConeTime += tocq(&coneTimer);
	return 0;
}
