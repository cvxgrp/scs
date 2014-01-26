#include "linAlg.h"
#include <math.h> 
/*
 * All basic linear operations are inlined (and further optimized) by the 
 * compiler. If compiling without optimization, causes code bloat.
 */

inline void _accumByAtrans(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y);
inline void _accumByA(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y);

// x = b*a
inline void setAsScaledArray(pfloat *x, const pfloat * a,const pfloat b,idxint len) {
  idxint i;
  for( i=0;i<len;++i ) x[i] = b*a[i];
}

// a*= b
inline void scaleArray(pfloat * a,const pfloat b,idxint len){
  idxint i;
  for( i=0;i<len;++i) a[i]*=b;
}

// x'*y
inline pfloat innerProd(const pfloat * x, const pfloat * y, idxint len){
  idxint i;
  pfloat ip = 0.0;
  for ( i=0;i<len;++i){
    ip += x[i]*y[i];
  }
  return ip;
}

// ||v||_2^2
inline pfloat calcNormSq(const pfloat * v,idxint len){
  idxint i;
  pfloat nmsq = 0.0;
  for ( i=0;i<len;++i){
    nmsq += v[i]*v[i];
  }
  return nmsq;
}

// ||v||_2
inline pfloat calcNorm(const pfloat * v,idxint len){
  return sqrt(calcNormSq(v, len));
}

inline pfloat calcNormInf(const pfloat *a, idxint l){
	pfloat tmp, max = 0.0;
	idxint i;
	for (i=0; i<l; ++i){
		tmp = fabs(a[i]);
		if(tmp > max) max = tmp;
	}
	return max;
}

// saxpy a += sc*b
inline void addScaledArray(pfloat * a, const pfloat * b, idxint n, const pfloat sc){
  idxint i;
  for (i=0;i<n;++i){
    a[i] += sc*b[i];
  }
}

inline pfloat calcNormDiff(const pfloat *a, const pfloat *b, idxint l) {
    pfloat nmDiff = 0.0, tmp;
    idxint i;
    for ( i=0; i<l; ++i){
        tmp = (a[i] - b[i]);
		nmDiff += tmp * tmp;
	}  
    return sqrt(nmDiff);
}

inline pfloat calcNormInfDiff(const pfloat *a, const pfloat *b, idxint l) {
	pfloat tmp, max = 0.0;
	idxint i;
	for ( i=0; i<l; ++i){
		tmp = fabs(a[i] - b[i]);
		if(tmp > max) max = tmp;
	}
	return max;
}

inline void accumByAtrans(Data * d, const pfloat *x, pfloat *y) 
{
	_accumByAtrans(d->n, d->Ax, d->Ai, d->Ap, x, y); 
}
inline void accumByA(Data * d, const pfloat *x, pfloat *y) 
{
	_accumByA(d->n, d->Ax, d->Ai, d->Ap, x, y);
}

inline void _accumByAtrans(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y) 
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

inline void _accumByA(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y) 
{
/*y  = A*x 
  A in column compressed format  
  this parallelizes over columns and uses
  pragma atomic to prevent concurrent writes to y 
 */
  idxint p, j;
  idxint c1, c2;
  pfloat xj;
//#pragma omp parallel for private(p,c1,c2,xj) 
  for (j = 0 ; j < n ; j++)
  {
      xj = x[j];
      c1 = Ap[j]; c2 = Ap[j+1];
      for (p = c1 ; p < c2 ; p++)        
      {
//#pragma omp atomic
          y [Ai[p]] += Ax [p] * xj ;
      }
  }
}
