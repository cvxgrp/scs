#ifndef SCS_MATRIX_H_GUARD
#define SCS_MATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

/** this struct defines the data matrices and is supplied in column compressed 
 * format: https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html
 */
typedef struct {
  /** matrix values, size: NNZ  */
  scs_float *x; 
  /** matrix row index, size: NNZ */
  scs_int *i;   
  /** matrix column pointer, size: n+1 */
  scs_int *p;   
  /** number of rows */
  scs_int m;
  /** number of columns */
  scs_int n; 
} ScsMatrix;

void SCS(_accum_by_atrans)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                           const scs_float *x, scs_float *y);
void SCS(_accum_by_a)(scs_int n, scs_float *Ax, scs_int *Ai, scs_int *Ap,
                      const scs_float *x, scs_float *y, scs_int skip_diag);
scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);

#ifdef __cplusplus
}
#endif
#endif
