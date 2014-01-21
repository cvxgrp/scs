#include "linAlg.h"
#include <math.h> 
/*
 * All basic linear operations are inlined (and further optimized) by the 
 * compiler. If compiling without optimization, causes code bloat.
 */

inline void _accumByAtrans(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y);
inline void _accumByA(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y);

// x = b*a
inline void setAsScaledArray(double *x, const double * a,const double b,int len) {
  int i;
  for( i=0;i<len;++i ) x[i] = b*a[i];
}

// a*= b
inline void scaleArray(double * a,const double b,int len){
  int i;
  for( i=0;i<len;++i) a[i]*=b;
}

// x'*y
inline double innerProd(const double * x, const double * y, int len){
  int i;
  double ip = 0.0;
  for ( i=0;i<len;++i){
    ip += x[i]*y[i];
  }
  return ip;
}

// ||v||_2^2
inline double calcNormSq(const double * v,int len){
  int i;
  double nmsq = 0.0;
  for ( i=0;i<len;++i){
    nmsq += v[i]*v[i];
  }
  return nmsq;
}

// ||v||_2
inline double calcNorm(const double * v,int len){
  return sqrt(calcNormSq(v, len));
}

inline double calcNormInf(const double *a, int l){
	double tmp, max = 0.0;
	int i;
	for (i=0; i<l; ++i){
		tmp = fabs(a[i]);
		if(tmp > max) max = tmp;
	}
	return max;
}

// saxpy a += sc*b
inline void addScaledArray(double * a, const double * b, int n, const double sc){
  int i;
  for (i=0;i<n;++i){
    a[i] += sc*b[i];
  }
}

inline double calcNormDiff(const double *a, const double *b, int l) {
    double nmDiff = 0.0, tmp;
    int i;
    for ( i=0; i<l; ++i){
        tmp = (a[i] - b[i]);
		nmDiff += tmp * tmp;
	}  
    return sqrt(nmDiff);
}

inline double calcNormInfDiff(const double *a, const double *b, int l) {
	double tmp, max = 0.0;
	int i;
	for ( i=0; i<l; ++i){
		tmp = fabs(a[i] - b[i]);
		if(tmp > max) max = tmp;
	}
	return max;
}

inline void accumByAtrans(Data * d, const double *x, double *y) 
{
	_accumByAtrans(d->n, d->Ax, d->Ai, d->Ap, x, y); 
}
inline void accumByA(Data * d, const double *x, double *y) 
{
	_accumByA(d->n, d->Ax, d->Ai, d->Ap, x, y);
}

inline void _accumByAtrans(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y) 
{
    /* y  = A'*x 
    A in column compressed format 
    parallelizes over columns (rows of A')
    */
    int p, j;
    int c1, c2; 
    double yj; 
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

inline void _accumByA(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y) 
{
/*y  = A*x 
  A in column compressed format  
  this parallelizes over columns and uses
  pragma atomic to prevent concurrent writes to y 
 */
  int p, j;
  int c1, c2;
  double xj;
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
