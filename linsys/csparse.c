/* Routines modified from CSparse, T. Davis et al */

#include "csparse.h"

csc *SCS(cs_spalloc)(scs_int m, scs_int n, scs_int nzmax, scs_int values,
                     scs_int triplet) {
  csc *A = (csc *)scs_calloc(1, sizeof(csc)); /* allocate the csc struct */
  if (!A) {
    return SCS_NULL;
  }         /* out of memory */
  A->m = m; /* define dimensions and nzmax */
  A->n = n;
  A->nzmax = nzmax = MAX(nzmax, 1);
  A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
  A->p = (scs_int *)scs_calloc((triplet ? nzmax : n + 1), sizeof(scs_int));
  A->i = (scs_int *)scs_calloc(nzmax, sizeof(scs_int));
  A->x = values ? (scs_float *)scs_calloc(nzmax, sizeof(scs_float)) : SCS_NULL;
  return (!A->p || !A->i || (values && !A->x)) ? SCS(cs_spfree)(A) : A;
}

csc *SCS(cs_done)(csc *C, void *w, void *x, scs_int ok) {
  scs_free(w); /* free workspace */
  scs_free(x);
  return ok ? C : SCS(cs_spfree)(C); /* return result if OK, else free it */
}

/* C = compressed-column form of a triplet matrix T */
csc *SCS(cs_compress)(const csc *T, scs_int *idx_mapping) {
  scs_int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  scs_float *Cx, *Tx;
  csc *C;
  m = T->m;
  n = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
  C = SCS(cs_spalloc)(m, n, nz, Tx != SCS_NULL, 0); /* allocate result */
  w = (scs_int *)scs_calloc(n, sizeof(scs_int));    /* get workspace */
  if (!C || !w) {
    return SCS(cs_done)(C, w, SCS_NULL, 0);
  } /* out of memory */
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (k = 0; k < nz; k++) {
    w[Tj[k]]++; /* column counts */
  }
  SCS(cumsum)(Cp, w, n); /* column pointers */
  for (k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
    if (idx_mapping) {
      idx_mapping[k] = p;
    } /* old to new indices */
    if (Cx) {
      Cx[p] = Tx[k];
    }
  }
  return SCS(cs_done)(C, w, SCS_NULL, 1); /* success; free w and return C */
}

scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n) {
  scs_int i, nz = 0;
  scs_float nz2 = 0;
  if (!p || !c) {
    return (-1);
  } /* check inputs */
  for (i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i]; /* also in scs_float to avoid scs_int overflow */
    c[i] = p[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p[n] = nz;
  return nz2; /* return sum (c [0..n-1]) */
}

csc *SCS(cs_spfree)(csc *A) {
  if (!A) {
    return SCS_NULL;
  } /* do nothing if A already SCS_NULL */
  scs_free(A->p);
  scs_free(A->i);
  scs_free(A->x);
  scs_free(A);
  return (csc *)SCS_NULL; /* free the csc struct and return SCS_NULL */
}
