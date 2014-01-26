#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD
#include "scs.h"
#include <math.h> 

void _accumByAtrans(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y);
void _accumByA(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y);
void setAsScaledArray(pfloat *x, const pfloat * a,const pfloat b,idxint len);
void scaleArray(pfloat * a,const pfloat b,idxint len);
pfloat innerProd(const pfloat * x, const pfloat * y, idxint len);
pfloat calcNormSq(const pfloat * v,idxint len);
pfloat calcNorm(const pfloat * v,idxint len);
pfloat calcNormInf(const pfloat *a, idxint l);
void addScaledArray(pfloat * a, const pfloat * b, idxint n, const pfloat sc);
pfloat calcNormDiff(const pfloat *a, const pfloat *b, idxint l);
pfloat calcNormInfDiff(const pfloat *a, const pfloat *b, idxint l);
void accumByAtrans(Data * d, const pfloat *x, pfloat *y);
void accumByA(Data * d, const pfloat *x, pfloat *y);
void _accumByAtrans(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y); 
void _accumByA(idxint n, pfloat * Ax, idxint * Ai, idxint * Ap, const pfloat *x, pfloat *y); 
#endif
