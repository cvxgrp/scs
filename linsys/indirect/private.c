#include "private.h"
#include "linAlg.h"

#define CG_BEST_TOL 1e-9
#define CG_EXPONENT 1.5
#define PRINT_INTERVAL 100

static void calcAx(Data * d, Priv * p, const pfloat * x, pfloat * y);
static idxint cgCustom(Data *d, Priv * p, const pfloat *s, pfloat * b, idxint max_its, pfloat tol);
static void transpose (Data * d, Priv * p);

static idxint totCgIts = 0;
static idxint lastNCgIts = 0;
static struct timeval tic_linsys_start;
static pfloat totalSolveTime;

void lTic(void) {
    gettimeofday(&tic_linsys_start, NULL);
}
pfloat lTocq(void) {
    struct timeval tic_timestop;
    gettimeofday(&tic_timestop, NULL);
    return tic_timestop.tv_sec*1e3 + tic_timestop.tv_usec/1e3 - tic_linsys_start.tv_sec*1e3 - tic_linsys_start.tv_usec/1e3;
}

char * getLinSysMethod() {
    char str[64];
    idxint len = sprintf(str,"sparse-indirect, CG tol ~ 1/iter^(%2.2f)", (pfloat) CG_EXPONENT);
    return strndup(str, len); 
}

char * getLinSysSummary(Info * info) {
    char str[64];
    idxint len = sprintf(str, "Avg num CG iterations: %2.2f, avg solve time %1.2es\n", (pfloat) totCgIts / (info->iter + 1), totalSolveTime / (info->iter + 1) / 1e3);
    totCgIts = 0;
    lastNCgIts = 0;
    totalSolveTime = 0;
    return strndup(str, len);
}

Priv * initPriv(Data * d) {
  Priv * p = scs_calloc(1,sizeof(Priv));
  p->p = scs_malloc((d->n)*sizeof(pfloat));
  p->r = scs_malloc((d->n)*sizeof(pfloat));
  p->Ap = scs_malloc((d->n)*sizeof(pfloat));
  p->tmp = scs_malloc((d->m)*sizeof(pfloat));

  p->Ati = scs_malloc((d->Ap[d->n])*sizeof(idxint));
  p->Atp = scs_malloc((d->m+1)*sizeof(idxint));
  p->Atx = scs_malloc((d->Ap[d->n])*sizeof(pfloat));
  transpose(d,p);
  totalSolveTime = 0;
  if (!p->p || !p->r || !p->Ap || !p->tmp || !p->Ati || !p->Atp || !p->Atx){
    freePriv(p);
    return NULL;
  }
  return p;
}

static void transpose (Data * d, Priv * p) {
  idxint * Ci = p->Ati;
  idxint * Cp = p->Atp;
  pfloat * Cx = p->Atx;
  idxint m = d->m;
  idxint n = d->n;

  idxint * Ap = d->Ap;
  idxint * Ai = d->Ai;
  pfloat * Ax = d->Ax;

  idxint i, j, q, *z;
  z = scs_calloc(m,sizeof(idxint));
  for (i = 0 ; i < Ap [n] ; i++) z [Ai [i]]++ ;          /* row counts */
  cs_cumsum (Cp, z, m) ;                                 /* row pointers */
  for (j = 0 ; j < n ; j++)
  {
    for (i = Ap [j] ; i < Ap [j+1] ; i++)
    {
      Ci [q = z [Ai [i]]++] = j ; /* place A(i,j) as entry C(j,i) */
      if (Cx) Cx [q] = Ax [i] ;
    }
  }
  scs_free(z);
}

void freePriv(Priv * p){
    if (p){
        if(p->p) scs_free(p->p);
        if(p->r) scs_free(p->r);
        if(p->Ap) scs_free(p->Ap);
        if(p->tmp) scs_free(p->tmp);
        if(p->Ati) scs_free(p->Ati);
        if(p->Atx) scs_free(p->Atx);
        if(p->Atp) scs_free(p->Atp);
        scs_free(p);
    }
}

void solveLinSys(Data *d, Priv * p, pfloat * b, const pfloat * s, idxint iter){
    idxint cgIts;
    pfloat cgTol = iter < 0 ? CG_BEST_TOL : calcNorm(b,d->n) / POWF(iter + 1, (pfloat) CG_EXPONENT);
	lTic();
    /* solves Mx = b, for x but stores result in b */
	/* s contains warm-start (if available) */
	accumByAtrans(d, p, &(b[d->n]), b);
   	/* solves (I+A'A)x = b, s warm start, solution stored in b */
	cgIts = cgCustom(d, p, s, b, d->n, cgTol);
    scaleArray(&(b[d->n]),-1,d->m);
	accumByA(d, p, b, &(b[d->n]));
	
    if(iter >= 0) {
        totCgIts += cgIts;
        lastNCgIts += cgIts;
        #ifdef EXTRAVERBOSE
        if (d->VERBOSE && (iter + 1) % PRINT_INTERVAL == 0) {
            scs_printf("\taverage CG iterations for last %i scs iters: %2.2f\n", (int) PRINT_INTERVAL, (pfloat) lastNCgIts / PRINT_INTERVAL);
            lastNCgIts = 0;
        }
        #endif
    }
    totalSolveTime += lTocq();
}

static idxint cgCustom(Data *d, Priv * pr, const pfloat * s, pfloat * b, idxint max_its, pfloat tol){
	/* solves (I+A'A)x = b */
	/* warm start cg with s */  
	idxint i = 0, n = d->n;
	pfloat rsold;
    pfloat *p = pr->p; /* cg direction */
	pfloat *Ap = pr->Ap; /* updated CG direction */
	pfloat *r = pr->r; /* cg residual */

	pfloat alpha, rsnew=0;
	if (s==NULL){
		memcpy(r,b,n*sizeof(pfloat));
		memset(b,0.0,n*sizeof(pfloat));
	}
	else{
		calcAx(d,pr,s,r);
		addScaledArray(r,b,n,-1);
		scaleArray(r,-1,n);
		memcpy(b,s,n*sizeof(pfloat));
	}
	memcpy(p,r,n*sizeof(pfloat));
	rsold=calcNorm(r,n);

	for (i=0; i < max_its; ++i){
		calcAx(d,pr,p,Ap);
		
		alpha=(rsold*rsold)/innerProd(p,Ap,n);
		addScaledArray(b,p,n,alpha);
		addScaledArray(r,Ap,n,-alpha);   

        rsnew=calcNorm(r,n);
        if (rsnew < tol){
            /*scs_printf("tol: %.4e, resid: %.4e, iters: %i\n", tol, rsnew, i+1); */
            return i+1;
        }   
        scaleArray(p,(rsnew*rsnew)/(rsold*rsold),n);
        addScaledArray(p,r,n,1);
        rsold=rsnew;
    }
    return i;
}

static void calcAx(Data * d, Priv * p, const pfloat * x, pfloat * y){
	pfloat * tmp = p->tmp;
	memset(tmp,0,d->m*sizeof(pfloat));
	accumByA(d,p,x,tmp);
	memset(y,0,d->n*sizeof(pfloat));
	accumByAtrans(d, p, tmp, y);
	addScaledArray(y,x,d->n,d->RHO_X);
}

void _accumByAtrans(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y) 
{
    /* y  = A'*x 
    A in column compressed format 
    parallelizes over columns (rows of A')
    */
    idxint p, j;
    idxint c1, c2; 
    pfloat yj; 
#pragma omp parallel for private(p,c1,c2,yj) 
    for (j = 0 ; j < n ; j++)
    {   
        yj = y[j];
        c1 = Ap[j]; c2 = Ap[j+1];
        for (p = c1 ; p < c2 ; p++)    
        {   
            yj += Ax[p] * x[ Ai[p] ] ; 
        }   
        y[j] = yj; 
    }   
}

void _accumByA(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y) 
{
/*y  = A*x 
  A in column compressed format  
  this parallelizes over columns and uses
  pragma atomic to prevent concurrent writes to y 
 */
  idxint p, j;
  idxint c1, c2;
  pfloat xj;
/*#pragma omp parallel for private(p,c1,c2,xj)  */
  for (j = 0 ; j < n ; j++)
  {
      xj = x[j];
      c1 = Ap[j]; c2 = Ap[j+1];
      for (p = c1 ; p < c2 ; p++)        
      {
/*#pragma omp atomic */
          y [Ai[p]] += Ax [p] * xj ;
      }
  }
}

void accumByAtrans(Data * d, Priv * p, const pfloat *x, pfloat *y) 
{
	_accumByAtrans(d->n, d->Ax, d->Ai, d->Ap, x, y); 
}
void accumByA(Data * d, Priv * p, const pfloat *x, pfloat *y) 
{
	_accumByAtrans(d->m, p->Atx, p->Ati, p->Atp, x, y);
}
