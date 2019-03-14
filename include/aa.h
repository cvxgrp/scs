#ifndef AA_H_GUARD
#define AA_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "glbopts.h"

typedef scs_float aa_float;
typedef scs_int aa_int;

typedef struct ACCEL_WORK AaWork;

/* Initialize Anderson Acceleration, allocates memory.
 *
 * Args:
 *  dim: the dimension of the variable for aa
 *  aa_mem: the memory (number of past iterations used) for aa
 *  type1: bool, if True use type 1 aa, otherwise use type 2
 *
 * Reurns:
 *  Pointer to aa workspace
 */
AaWork *aa_init(aa_int dim, aa_int aa_mem, aa_int type1);

/* Apply Anderson Acceleration.
 *
 * Args:
 *  f: output of map at current iteration, overwritten with aa output at end.
 *  x: input to map at current iteration
 *  a: aa workspace from aa_init
 *
 * Returns:
 *  int, a value of 0 is success, <0 is failure at which point f is unchanged
 */
aa_int aa_apply(aa_float *f, const aa_float *x, AaWork *a);

/* Finish Anderson Acceleration, clears memory.
 *
 * Args:
 *  a: aa workspace from aa_init.
 */
void aa_finish(AaWork *a);

#define MAX_AA_NRM (1e4)

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifdef __cplusplus
}
#endif
#endif
