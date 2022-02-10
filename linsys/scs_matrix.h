#ifndef SCS_MATRIX_H_GUARD
#define SCS_MATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"
#include "scs_work.h"

/* Normalization routines, used if d->NORMALIZE is true */
/* normalizes A matrix, sets scal->E and scal->D diagonal scaling matrices,
 * A -> D*A*E. D and E must be all positive entries, D must satisfy cone
 * boundaries */
ScsScaling *SCS(normalize_a_p)(ScsMatrix *P, ScsMatrix *A, ScsConeWork *cone);

/* unnormalizes A matrix, unnormalizes by w->D and w->E */
/* void SCS(un_normalize_a_p)(ScsMatrix *A, ScsMatrix *P, const ScsScaling
 * *scal);
 */

/* to free the memory allocated in a ScsMatrix (called on A and P at finish) */
void SCS(free_scs_matrix)(ScsMatrix *A);

/* copies A (instead of in-place normalization), returns 0 for failure,
 * allocates memory for dstp	*/
scs_int SCS(copy_matrix)(ScsMatrix **dstp, const ScsMatrix *src);

scs_float SCS(cumsum)(scs_int *p, scs_int *c, scs_int n);

/**
 * Validate the linear system inputs, returns < 0 if not valid inputs.
 *
 *  @param  A    A data matrix
 *  @param  P    P data matrix
 *  @return status < 0 indicates failure
 */
scs_int SCS(validate_lin_sys)(const ScsMatrix *A, const ScsMatrix *P);

/**
 * Forms y += A^T * x
 *
 *  @param  A    A data matrix
 *  @param  x    Input
 *  @param  y    Output
 */
void SCS(accum_by_atrans)(const ScsMatrix *A, const scs_float *x, scs_float *y);

/**
 * Forms y += A * x
 *
 *  @param  A           Data matrix
 *  @param  x           Input
 *  @param  y           Output
 */
void SCS(accum_by_a)(const ScsMatrix *A, const scs_float *x, scs_float *y);

/**
 * Forms y += P * x
 *
 *  @param  P    P data matrix
 *  @param  x    Input
 *  @param  y    Output
 */
void SCS(accum_by_p)(const ScsMatrix *P, const scs_float *x, scs_float *y);

#ifdef __cplusplus
}
#endif
#endif
