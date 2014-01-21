#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD
#include "scs.h"
#include <math.h> 

void _accumByAtrans(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y);
void _accumByA(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y);
void setAsScaledArray(double *x, const double * a,const double b,int len);
void scaleArray(double * a,const double b,int len);
double innerProd(const double * x, const double * y, int len);
double calcNormSq(const double * v,int len);
double calcNorm(const double * v,int len);
double calcNormInf(const double *a, int l);
void addScaledArray(double * a, const double * b, int n, const double sc);
double calcNormDiff(const double *a, const double *b, int l);
double calcNormInfDiff(const double *a, const double *b, int l);
void accumByAtrans(Data * d, const double *x, double *y);
void accumByA(Data * d, const double *x, double *y);
void _accumByAtrans(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y); 
void _accumByA(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y); 
#endif
