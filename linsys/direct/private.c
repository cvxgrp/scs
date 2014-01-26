#include "private.h"

// forward declare
idxint LDLInit(cs * A, idxint P[], pfloat **info);
idxint LDLFactor(cs * A, idxint P[], idxint Pinv[], cs ** L, pfloat **D);
void LDLSolve(pfloat *x, pfloat b[], cs * L, pfloat D[], idxint P[], pfloat * bp);
idxint factorize(Data * d,Work * w);

char * getLinSysSummary(Data * d, Info * info) {
    return NULL;
}

void freePriv(Work * w){
	if(w) {
        if(w->p) {
            if (w->p->L) cs_spfree(w->p->L);
            if (w->p->P) scs_free(w->p->P);
            if (w->p->D) scs_free(w->p->D);
            if (w->p->bp) scs_free(w->p->bp);
            scs_free(w->p);
	    }
    }
}

void solveLinSys(Data *d, Work * w, pfloat * b, const pfloat * s, idxint iter){
  // returns solution to linear system
  // Ax = b with solution stored in b
  LDLSolve(b, b, w->p->L, w->p->D, w->p->P, w->p->bp);
}

idxint privateInitWork(Data * d, Work * w){ 
	w->method = strdup("sparse-direct");
	idxint n_plus_m = d->n + d->m;
	w->p = scs_malloc(sizeof(Priv));
	w->p->P = scs_malloc(sizeof(idxint)*n_plus_m);
	w->p->L = scs_malloc(sizeof (cs));
	w->p->bp = scs_malloc(n_plus_m * sizeof(pfloat));
	w->p->L->m = n_plus_m;
	w->p->L->n = n_plus_m;
	w->p->L->nz = -1; 
	return(factorize(d,w));
}

cs * formKKT(Data * d, Work * w){
	/* ONLY UPPER TRIANGULAR PART IS STUFFED
	 * forms column compressed KKT matrix
	 * assumes column compressed form A matrix
	 *
	 * forms upper triangular part of [I A'; A -I]
	 */
	idxint j, k, kk;
	/* I at top left */
	const idxint Anz = d->Ap[d->n];
	const idxint Knzmax = d->n + d->m + Anz;
	cs * K = cs_spalloc(d->m + d->n, d->m + d->n, Knzmax, 1, 1);
	if (!K)
    {
        return NULL;
    }
    kk = 0;
	for (k = 0; k < d->n; k++){
		K->i[kk] = k;
		K->p[kk] = k;
		K->x[kk] = d->RHO_X;
		kk++;
	}
	/* A^T at top right : CCS: */
	for (j = 0; j < d->n; j++) {                 
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			K->p[kk] = d->Ai[k] + d->n;
			K->i[kk] = j;
			K->x[kk] = d->Ax[k];
			kk++;
		}   
	}
	/* -I at bottom right */
	for (k = 0; k < d->m; k++){
		K->i[kk] = k + d->n;
		K->p[kk] = k + d->n;
		K->x[kk] = -1;
		kk++;
	}
	// assert kk == Knzmax
	K->nz = Knzmax;
	cs * K_cs = cs_compress(K);
	cs_spfree(K);
	return(K_cs);
}


idxint factorize(Data * d,Work * w){
	//tic();
	cs * K = formKKT(d,w);
	if (!K){
        return -7; //arbitrary int
    }
    //if(d->VERBOSE) scs_printf("KKT matrix factorization info:\n");
	pfloat *info;
	idxint amd_status = LDLInit(K, w->p->P, &info);
	if (amd_status < 0) return(amd_status);
	/*
  if(d->VERBOSE) {
#ifdef DLONG
		amd_l_info(info);
#else
		amd_info(info);
#endif
	}
  */
	idxint * Pinv = cs_pinv(w->p->P, w->l-1);
	cs * C = cs_symperm(K, Pinv, 1); 
	idxint ldl_status = LDLFactor(C, NULL, NULL, &w->p->L, &w->p->D);
	//if(d->VERBOSE) scs_printf("KKT matrix factorization took %4.8f ms\n",tocq());
	cs_spfree(C);cs_spfree(K);scs_free(Pinv);scs_free(info);
	return(ldl_status);
}

idxint LDLInit(cs * A, idxint P[], pfloat **info) {
	*info  = (pfloat *) scs_malloc(AMD_INFO * sizeof(pfloat));
#ifdef DLONG
	return(amd_l_order(A->n, A->p, A->i, P, (pfloat *) NULL, *info));
#else
	return(amd_order(A->n, A->p, A->i, P, (pfloat *) NULL, *info));
#endif
}

idxint LDLFactor(cs * A, idxint P[], idxint Pinv[], cs **L , pfloat **D) 
{
	idxint n = A->n;
	(*L)->p = (idxint *) scs_malloc((1 + n) * sizeof(idxint));
	
	//idxint Parent[n], Lnz[n], Flag[n], Pattern[n];
	//pfloat Y[n];

	idxint * Parent = scs_malloc(n*sizeof(idxint));
	idxint * Lnz = scs_malloc(n*sizeof(idxint));
	idxint * Flag = scs_malloc(n*sizeof(idxint));
	idxint * Pattern = scs_malloc(n*sizeof(idxint));
	pfloat * Y = scs_malloc(n*sizeof(pfloat));

	LDL_symbolic(n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);

	(*L)->nzmax = *((*L)->p + n);
	(*L)->x = (pfloat *) scs_malloc((*L)->nzmax * sizeof(pfloat));
	(*L)->i =    (idxint *) scs_malloc((*L)->nzmax * sizeof(idxint));
	*D  = (pfloat *) scs_malloc(n * sizeof(pfloat));
	
	if(!(*D) || !(*L)->i || !(*L)->x || !Y || !Pattern || !Flag || !Lnz || !Parent) return -1;
	
	idxint kk = LDL_numeric(n, A->p, A->i, A->x, (*L)->p, Parent, Lnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);

	scs_free(Parent);
	scs_free(Lnz);
	scs_free(Flag);
	scs_free(Pattern);
	scs_free(Y);
	return(n - kk);
}

void LDLSolve(pfloat *x, pfloat b[], cs * L, pfloat D[], idxint P[], pfloat * bp)
{
  // solves PLDL'P' x = b for x
  idxint n = L->n;
  if (P == NULL) {
    if (x != b) // if they're different addresses
      memcpy(x,b, n*sizeof(pfloat)); 
    LDL_lsolve(n, x, L->p, L->i, L->x);
    LDL_dsolve(n, x, D); 
    LDL_ltsolve(n, x, L->p, L->i, L->x);
  } else {
    LDL_perm(n, bp, b, P); 
    LDL_lsolve(n, bp, L->p, L->i, L->x);
    LDL_dsolve(n, bp, D); 
    LDL_ltsolve(n, bp, L->p, L->i, L->x);
    LDL_permt(n, x, bp, P); 
	}   
}
