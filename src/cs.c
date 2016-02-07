#include "cs.h"
#include "scs.h"
/* NB: this is a subset of the routines in the CSPARSE package by
 Tim Davis et. al., for the full package please visit
 http://www.cise.ufl.edu/research/sparse/CSparse/ */

/* wrapper for malloc */
static void *cs_malloc(scs_int n, scs_int size) {
    return (scs_malloc(n * size));
}

/* wrapper for calloc */
static void *cs_calloc(scs_int n, scs_int size) {
    return (scs_calloc(n, size));
}

/* wrapper for free */
static void *cs_free(void *p) {
    if (p)
        scs_free(p);   /* free p if it is not already SCS_NULL */
    return (SCS_NULL); /* return SCS_NULL to simplify the use of cs_free */
}

/* C = compressed-column form of a triplet matrix T */
cs *cs_compress(const cs *T) {
    scs_int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
    scs_float *Cx, *Tx;
    cs *C;
    m = T->m;
    n = T->n;
    Ti = T->i;
    Tj = T->p;
    Tx = T->x;
    nz = T->nz;
    C = cs_spalloc(m, n, nz, Tx != SCS_NULL, 0); /* allocate result */
    w = cs_calloc(n, sizeof(scs_int));           /* get workspace */
    if (!C || !w)
        return (cs_done(C, w, SCS_NULL, 0)); /* out of memory */
    Cp = C->p;
    Ci = C->i;
    Cx = C->x;
    for (k = 0; k < nz; k++)
        w[Tj[k]]++;      /* column counts */
    cs_cumsum(Cp, w, n); /* column pointers */
    for (k = 0; k < nz; k++) {
        Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
        if (Cx)
            Cx[p] = Tx[k];
    }
    return (cs_done(C, w, SCS_NULL, 1)); /* success; free w and return C */
}

cs *cs_done(cs *C, void *w, void *x, scs_int ok) {
    cs_free(w); /* free workspace */
    cs_free(x);
    return (ok ? C : cs_spfree(C)); /* return result if OK, else free it */
}

cs *cs_spalloc(scs_int m, scs_int n, scs_int nzmax, scs_int values,
               scs_int triplet) {
    cs *A = cs_calloc(1, sizeof(cs)); /* allocate the cs struct */
    if (!A)
        return (SCS_NULL); /* out of memory */
    A->m = m;              /* define dimensions and nzmax */
    A->n = n;
    A->nzmax = nzmax = MAX(nzmax, 1);
    A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
    A->p = cs_malloc(triplet ? nzmax : n + 1, sizeof(scs_int));
    A->i = cs_malloc(nzmax, sizeof(scs_int));
    A->x = values ? cs_malloc(nzmax, sizeof(scs_float)) : SCS_NULL;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree(A) : A);
}

cs *cs_spfree(cs *A) {
    if (!A)
        return (SCS_NULL); /* do nothing if A already SCS_NULL */
    cs_free(A->p);
    cs_free(A->i);
    cs_free(A->x);
    return ((cs *)cs_free(A)); /* free the cs struct and return SCS_NULL */
}

scs_float cs_cumsum(scs_int *p, scs_int *c, scs_int n) {
    scs_int i, nz = 0;
    scs_float nz2 = 0;
    if (!p || !c)
        return (-1); /* check inputs */
    for (i = 0; i < n; i++) {
        p[i] = nz;
        nz += c[i];
        nz2 += c[i]; /* also in scs_float to avoid scs_int overflow */
        c[i] = p[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p[n] = nz;
    return (nz2); /* return sum (c [0..n-1]) */
}

scs_int *cs_pinv(scs_int const *p, scs_int n) {
    scs_int k, *pinv;
    if (!p)
        return (SCS_NULL);                /* p = SCS_NULL denotes identity */
    pinv = cs_malloc(n, sizeof(scs_int)); /* allocate result */
    if (!pinv)
        return (SCS_NULL); /* out of memory */
    for (k = 0; k < n; k++)
        pinv[p[k]] = k; /* invert the permutation */
    return (pinv);      /* return result */
}

cs *cs_symperm(const cs *A, const scs_int *pinv, scs_int values) {
    scs_int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
    scs_float *Cx, *Ax;
    cs *C;
    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    C = cs_spalloc(n, n, Ap[n], values && (Ax != SCS_NULL),
                   0);                 /* alloc result*/
    w = cs_calloc(n, sizeof(scs_int)); /* get workspace */
    if (!C || !w)
        return (cs_done(C, w, SCS_NULL, 0)); /* out of memory */
    Cp = C->p;
    Ci = C->i;
    Cx = C->x;
    for (j = 0; j < n; j++) /* count entries in each column of C */
    {
        j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            if (i > j)
                continue;            /* skip lower triangular part of A */
            i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
            w[MAX(i2, j2)]++;        /* column count of C */
        }
    }
    cs_cumsum(Cp, w, n); /* compute column pointers of C */
    for (j = 0; j < n; j++) {
        j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            if (i > j)
                continue;            /* skip lower triangular part of A*/
            i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
            Ci[q = w[MAX(i2, j2)]++] = MIN(i2, j2);
            if (Cx)
                Cx[q] = Ax[p];
        }
    }
    return (cs_done(C, w, SCS_NULL, 1)); /* success; free workspace, return C */
}
