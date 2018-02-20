#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

/* YOUR LINEAR SYSTEM SOLVER MUST IMPLEMENT THESE METHODS AND SCS_LIN_SYS_WORK
 * STRUCT */

/* initialize ScsLinSysWork structure and perform any necessary preprocessing */
ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A,
                                      const ScsSettings *stgs);
/* solves [d->RHO_X * I  A' ; A  -I] x = b for x, stores result in b, s contains
 * warm-start, iter is current scs iteration count */
scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsSettings *stgs,
                           ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_int iter);
/* frees ScsLinSysWork structure and allocated memory in ScsLinSysWork */
void SCS(free_lin_sys_work)(ScsLinSysWork *p);

/* forms y += A'*x */
void SCS(accum_by_atrans)(const ScsMatrix *A, ScsLinSysWork *p,
                          const scs_float *x, scs_float *y);
/* forms y += A*x */
void SCS(accum_by_a)(const ScsMatrix *A, ScsLinSysWork *p, const scs_float *x,
                     scs_float *y);

/* returns negative num if input data is invalid */
scs_int SCS(validate_lin_sys)(const ScsMatrix *A);

/* returns string describing method, can return null, if not null free will be
 * called on output */
char *SCS(get_lin_sys_method)(const ScsMatrix *A, const ScsSettings *stgs);
/* returns string containing summary information about linear system solves, can
 * return null, if not null free will be called on output */
char *SCS(get_lin_sys_summary)(ScsLinSysWork *p, const ScsInfo *info);

/* Normalization routines, used if d->NORMALIZE is true */
/* normalizes A matrix, sets w->E and w->D diagonal scaling matrices, Anew =
 * d->SCALE * (D^-1)*A*(E^-1) (different to paper which is D*A*E)
 * D and E must be all positive entries, D must satisfy cone boundaries
 * must set (w->mean_norm_row_a = mean of norms of rows of normalized A) THEN
 * scale
 * resulting A by d->SCALE */
void SCS(normalize_a)(ScsMatrix *A, const ScsSettings *stgs, const ScsCone *k,
                      ScsScaling *scal);
/* unnormalizes A matrix, unnormalizes by w->D and w->E and d->SCALE */
void SCS(un_normalize_a)(ScsMatrix *A, const ScsSettings *stgs,
                         const ScsScaling *scal);
/* to free the memory allocated in ScsMatrix */
void SCS(free_a_matrix)(ScsMatrix *A);

/* copies A (instead of in-place normalization), returns 0 for failure,
 * allocates memory for dstp	*/
scs_int SCS(copy_a_matrix)(ScsMatrix **dstp, const ScsMatrix *src);

#ifdef __cplusplus
}
#endif

#endif
