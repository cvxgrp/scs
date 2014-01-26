#include "private.h"
#include "linAlg.h"

#define CG_BEST_TOL 1e-9
#define CG_EXPONENT 1.5
#define CG_VERBOSE 0
#define PRINT_INTERVAL 100

static inline void calcAx(Data * d, Work * w, const pfloat * x, pfloat * y);
static idxint cgCustom(Data *d, Work *w, const pfloat *s, pfloat * b, idxint max_its, pfloat tol);
static inline void CGaccumByA(Data * d, Work * w, const pfloat *x, pfloat *y);
static inline void CGaccumByAtrans(Data *d, Work * w, const pfloat *x, pfloat *y);
static inline void transpose (Data * d, Work * w);

static idxint totCgIts = 0;
static idxint lastNCgIts = 0;

char * getLinSysSummary(Data * d, Info * info) {
    char str[80];
    idxint len = sprintf(str, "Average CG iterations: %2.2f\n", (float) totCgIts / info->iter );
    return strndup(str, len);
}

idxint privateInitWork(Data * d, Work * w){
  char str[80];
  idxint len = sprintf(str,"sparse-indirect, CG tol ~ 1/iter^(%2.2f)", (float) CG_EXPONENT);
  w->method = strndup(str, len);
  w->p = scs_malloc(sizeof(Priv));
  w->p->p = scs_malloc((d->n)*sizeof(pfloat));
  w->p->r = scs_malloc((d->n)*sizeof(pfloat));
  w->p->Ap = scs_malloc((d->n)*sizeof(pfloat));
  w->p->tmp = scs_malloc((d->m)*sizeof(pfloat));

  w->p->Ati = scs_malloc((d->Ap[d->n])*sizeof(idxint));
  w->p->Atp = scs_malloc((d->m+1)*sizeof(idxint));
  w->p->Atx = scs_malloc((d->Ap[d->n])*sizeof(pfloat));
  transpose(d,w);
  return(0);
}

static inline void transpose (Data * d, Work * w)
{
  idxint * Ci = w->p->Ati;
  idxint * Cp = w->p->Atp;
  pfloat * Cx = w->p->Atx;
  idxint m = d->m;
  idxint n = d->n;

  idxint * Ap = d->Ap;
  idxint * Ai = d->Ai;
  pfloat * Ax = d->Ax;

  idxint p, j, q, *z;
  z = scs_calloc(m,sizeof(idxint));
  for (p = 0 ; p < Ap [n] ; p++) z [Ai [p]]++ ;          /* row counts */
  cs_cumsum (Cp, z, m) ;                                 /* row pointers */
  for (j = 0 ; j < n ; j++)
  {
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      Ci [q = z [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
      if (Cx) Cx [q] = Ax [p] ;
    }
  }
  scs_free(z);
}

void freePriv(Work * w){
    if (w) {
        if (w->p){
            if(w->p->p) scs_free(w->p->p);
            if(w->p->r) scs_free(w->p->r);
            if(w->p->Ap) scs_free(w->p->Ap);
            if(w->p->tmp) scs_free(w->p->tmp);
            if(w->p->Ati) scs_free(w->p->Ati);
            if(w->p->Atx) scs_free(w->p->Atx);
            if(w->p->Atp) scs_free(w->p->Atp);
            scs_free(w->p);
        }
    }
}

void solveLinSys(Data *d, Work * w, pfloat * b, const pfloat * s, idxint iter){
    pfloat cgTol = iter < 0 ? CG_BEST_TOL : calcNorm(b,d->n) / pow(iter + 1, (float) CG_EXPONENT);
	// solves Mx = b, for x but stores result in b
	// s contains warm-start (if available)
	CGaccumByAtrans(d,w, &(b[d->n]), b);
   	// solves (I+A'A)x = b, s warm start, solution stored in b
	idxint cgIts = cgCustom(d, w, s, b, d->n, cgTol);
    scaleArray(&(b[d->n]),-1,d->m);
	CGaccumByA(d, w, b, &(b[d->n]));
	
    if(iter >= 0) {
        totCgIts += cgIts;
        lastNCgIts += cgIts;
        if (CG_VERBOSE && d->VERBOSE && (iter + 1) % PRINT_INTERVAL == 0) {
            scs_printf("\taverage CG iterations for last %i scs iters: %2.2f\n", (int) PRINT_INTERVAL, (pfloat) lastNCgIts / PRINT_INTERVAL);
            lastNCgIts = 0;
        }
    }
}

static idxint cgCustom(Data *d, Work *w, const pfloat * s, pfloat * b, idxint max_its, pfloat tol){
	/* solves (I+A'A)x = b */
	/* warm start cg with s */  
	idxint i = 0, n = d->n;
	pfloat *p = w->p->p; // cg direction
	pfloat *Ap = w->p->Ap; // updated CG direction
	pfloat *r = w->p->r; // cg residual

	pfloat alpha, rsnew=0;
	if (s==NULL){
		memcpy(r,b,n*sizeof(pfloat));
		memset(b,0.0,n*sizeof(pfloat));
	}
	else{
		calcAx(d,w,s,r);
		addScaledArray(r,b,n,-1);
		scaleArray(r,-1,n);
		memcpy(b,s,n*sizeof(pfloat));
	}
	memcpy(p,r,n*sizeof(pfloat));
	pfloat rsold=calcNorm(r,n);

	for (i=0; i < max_its; ++i){
		calcAx(d,w,p,Ap);
		
		alpha=(rsold*rsold)/innerProd(p,Ap,n);
		addScaledArray(b,p,n,alpha);
		addScaledArray(r,Ap,n,-alpha);   

        rsnew=calcNorm(r,n);
        if (rsnew < tol){
            //scs_printf("tol: %.4e, resid: %.4e, iters: %i\n", tol, rsnew, i+1);
            return i+1;
        }   
        scaleArray(p,(rsnew*rsnew)/(rsold*rsold),n);
        addScaledArray(p,r,n,1);
        rsold=rsnew;
    }
    return i;
}

static inline void calcAx(Data * d, Work * w, const pfloat * x, pfloat * y){
	pfloat * tmp = w->p->tmp;
	memset(tmp,0,d->m*sizeof(pfloat));
	CGaccumByA(d,w,x,tmp);
	memset(y,0,d->n*sizeof(pfloat));
	CGaccumByAtrans(d,w,tmp,y);
	addScaledArray(y,x,d->n,d->RHO_X);
}

static inline void CGaccumByA(Data * d, Work * w, const pfloat *x, pfloat *y)
{
  _accumByAtrans(d->m,w->p->Atx,w->p->Ati,w->p->Atp,x,y);
}
static inline void CGaccumByAtrans(Data *d, Work * w, const pfloat *x, pfloat *y)
{
  _accumByAtrans(d->n,d->Ax,d->Ai,d->Ap,x,y);
}
