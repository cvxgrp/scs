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
 *  @param  rho_x      rho_x float
 *  @return            Linear system solver workspace
 *
 */
ScsLinSysWork *SCS(init_lin_sys_work)(const ScsMatrix *A, const ScsMatrix *P,
                                      scs_float *rho_y_vec, scs_float rho_x);

/**
 * Frees ScsLinSysWork structure and allocated memory in ScsLinSysWork
 *
 *  @param  w    Linear system private workspace
 */
void SCS(free_lin_sys_work)(ScsLinSysWork *w);

/**
 * Solves the linear system required by SCS at each iteration:
 * \f[
 *    \begin{bmatrix}
 *    (\rho_x I + P) & A^\top \\
 *     A   &  -\mathrm{diag}(\rho_y) \\
 *    \end{bmatrix} x = b
 *  \f]
 *
 *  for x, result stored result in b
 *
 *  @param  w    Linear system private workspace
 *  @param  b    Right hand side, should contain solution at the end
 *  @param  s    Contains optional warm-start
 *  @param  tol  Tolerance required for the system solve
 *  @return status < 0 indicates failure
 *
 */
scs_int SCS(solve_lin_sys)(ScsLinSysWork *w, scs_float *b, const scs_float *s,
                           scs_float tol);
/**
 *  Update the linsys workspace when new rho_y_vec is updated
 *
 *  @param  w          Linear system private workspace
 *  @param  rho_y_vec  R_y diagonal entries
 *
 */
void SCS(update_lin_sys_rho_y_vec)(ScsLinSysWork *w, scs_float *rho_y_vec);

/**
 * Name of the linear solver. Can return null.
 *
 * @return name of method
 */
const char *SCS(get_lin_sys_method)(void);

#ifdef __cplusplus
}
#endif

#endif
