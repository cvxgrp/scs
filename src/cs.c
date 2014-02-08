#include "cs.h"
#include "scs.h"
/* NB: this is a subset of the routines in the CSPARSE package by
 Tim Davis et. al., for the full package please visit
 http://www.cise.ufl.edu/research/sparse/CSparse/ */

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))

/* wrapper for malloc */
static void *cs_malloc(idxint n, idxint size) {
	return (scs_malloc(n * size));
}

/* wrapper for calloc */
static void *cs_calloc(idxint n, idxint size) {
	return (scs_calloc(n, size));
}

/* wrapper for free */
static void *cs_free(void *p) {
	if (p) scs_free(p); /* free p if it is not already NULL */
	return (NULL); /* return NULL to simplify the use of cs_free */
}

/* C = compressed-column form of a triplet matrix T */
cs *cs_compress(const cs *T) {
	idxint m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
	pfloat *Cx, *Tx;
	cs *C;
	m = T->m;
	n = T->n;
	Ti = T->i;
	Tj = T->p;
	Tx = T->x;
	nz = T->nz;
	C = cs_spalloc(m, n, nz, Tx != NULL, 0); /* allocate result */
	w = cs_calloc(n, sizeof(idxint)); /* get workspace */
	if (!C || !w)
		return (cs_done(C, w, NULL, 0)); /* out of memory */
	Cp = C->p;
	Ci = C->i;
	Cx = C->x;
	for (k = 0; k < nz; k++)
		w[Tj[k]]++; /* column counts */
	cs_cumsum(Cp, w, n); /* column pointers */
	for (k = 0; k < nz; k++) {
		Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
		if (Cx)
			Cx[p] = Tx[k];
	}
	return (cs_done(C, w, NULL, 1)); /* success; free w and return C */
}

cs *cs_done(cs *C, void *w, void *x, idxint ok) {
	cs_free(w); /* free workspace */
	cs_free(x);
	return (ok ? C : cs_spfree(C)); /* return result if OK, else free it */
}

cs *cs_spalloc(idxint m, idxint n, idxint nzmax, idxint values, idxint triplet) {
	cs *A = cs_calloc(1, sizeof(cs)); /* allocate the cs struct */
	if (!A)
		return (NULL); /* out of memory */
	A->m = m; /* define dimensions and nzmax */
	A->n = n;
	A->nzmax = nzmax = CS_MAX(nzmax, 1);
	A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
	A->p = cs_malloc(triplet ? nzmax : n + 1, sizeof(idxint));
	A->i = cs_malloc(nzmax, sizeof(idxint));
	A->x = values ? cs_malloc(nzmax, sizeof(pfloat)) : NULL;
	return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree(A) : A);
}

cs *cs_spfree(cs *A) {
	if (!A)
		return (NULL); /* do nothing if A already NULL */
	cs_free(A->p);
	cs_free(A->i);
	cs_free(A->x);
	return ((cs *) cs_free(A)); /* free the cs struct and return NULL */
}

pfloat cs_cumsum(idxint *p, idxint *c, idxint n) {
	idxint i, nz = 0;
	pfloat nz2 = 0;
	if (!p || !c)
		return (-1); /* check inputs */
	for (i = 0; i < n; i++) {
		p[i] = nz;
		nz += c[i];
		nz2 += c[i]; /* also in pfloat to avoid idxint overflow */
		c[i] = p[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
	}
	p[n] = nz;
	return (nz2); /* return sum (c [0..n-1]) */
}

idxint *cs_pinv(idxint const *p, idxint n) {
	idxint k, *pinv;
	if (!p)
		return (NULL); /* p = NULL denotes identity */
	pinv = cs_malloc(n, sizeof(idxint)); /* allocate result */
	if (!pinv)
		return (NULL); /* out of memory */
	for (k = 0; k < n; k++)
		pinv[p[k]] = k;/* invert the permutation */
	return (pinv); /* return result */
}

cs *cs_symperm(const cs *A, const idxint *pinv, idxint values) {
	idxint i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
	pfloat *Cx, *Ax;
	cs *C;
	n = A->n;
	Ap = A->p;
	Ai = A->i;
	Ax = A->x;
	C = cs_spalloc(n, n, Ap[n], values && (Ax != NULL), 0); /* alloc result*/
	w = cs_calloc(n, sizeof(idxint)); /* get workspace */
	if (!C || !w)
		return (cs_done(C, w, NULL, 0)); /* out of memory */
	Cp = C->p;
	Ci = C->i;
	Cx = C->x;
	for (j = 0; j < n; j++) /* count entries in each column of C */
	{
		j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */
		for (p = Ap[j]; p < Ap[j + 1]; p++) {
			i = Ai[p];
			if (i > j)
				continue; /* skip lower triangular part of A */
			i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
			w[CS_MAX(i2, j2)]++; /* column count of C */
		}
	}
	cs_cumsum(Cp, w, n); /* compute column pointers of C */
	for (j = 0; j < n; j++) {
		j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */
		for (p = Ap[j]; p < Ap[j + 1]; p++) {
			i = Ai[p];
			if (i > j)
				continue; /* skip lower triangular part of A*/
			i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
			Ci[q = w[CS_MAX(i2, j2)]++] = CS_MIN(i2, j2);
			if (Cx)
				Cx[q] = Ax[p];
		}
	}
	return (cs_done(C, w, NULL, 1)); /* success; free workspace, return C */
}
