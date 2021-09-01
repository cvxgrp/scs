#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

/* This is the API that any new linear system solver must implement */

/**
 * Initialize ScsLinSysWork structure and perform any necessary preprocessing
 *
 *  @param  A          A data matrix
 *  @param  P          P data matrix
 *  @param  rho_y_vec  R_y diagonal entries
 *  @param  rho_x SCS  rho_x float
 *  @return            Linear system solver workspace
 *
 */
ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A, const ScsMatrix *P,
                                      scs_float *rho_y_vec, scs_float rho_x);

/**
 * Solves the linear system required by SCS at each iteration:
 *
 *    [(rho_x * I + P)  A' ] x = b
 *    [ A              -R_y]
 *
 *  for x, result stored result in b
 *
 *  @param  A    A data matrix
 *  @param  P    P data matrix
 *  @param  p    Linear system private workspace
 *  @param  b    Right hand side, should contain solution at the end
 *  @param  s    Contains optional warm-start
 *  @param  tol  Tolerance required for the system solve
 *  @return status < 0 indicates failure
 *
 */
scs_int SCS(solve_lin_sys)(const ScsMatrix *A, const ScsMatrix *P,
                           ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol);

/**
 * Frees ScsLinSysWork structure and allocated memory in ScsLinSysWork
 *
 *  @param  p    Linear system private workspace
 */
void SCS(free_lin_sys_work)(ScsLinSysWork *p);


/**
 *  Update the linsys workspace when new rho_y_vec is updated 
 *
 *  @param  A          A data matrix
 *  @param  P          P data matrix
 *  @param  p          Linear system private workspace
 *  @param  rho_y_vec  R_y diagonal entries
 *
 */
void SCS(update_linsys_rho_y_vec)(const ScsMatrix *A, const ScsMatrix *P,
                                  ScsLinSysWork *p, scs_float *rho_y_vec);

/**
 * Name of the linear solver. Can return null. If not null free will be called
 * on the string.
 *
 * @return name of method
 */
char *SCS(get_lin_sys_method)(void);

#ifdef __cplusplus
}
#endif

#endif
