/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* Copyright (c) Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved.  See LDL/Doc/License.txt for the License.
 */

#include "SuiteSparse_config.h"
#define LDL_int scs_int
#define LDL_ID "%d"

#define LDL_symbolic ldl_symbolic
#define LDL_numeric ldl_numeric
#define LDL_lsolve ldl_lsolve
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve
#define LDL_perm ldl_perm
#define LDL_permt ldl_permt
#define LDL_valid_perm ldl_valid_perm
#define LDL_valid_matrix ldl_valid_matrix

void ldl_symbolic (scs_int n, scs_int Ap [ ], scs_int Ai [ ], scs_int Lp [ ],
    scs_int Parent [ ], scs_int Lnz [ ], scs_int Flag [ ], scs_int P [ ], scs_int Pinv [ ]) ;

scs_int ldl_numeric (scs_int n, scs_int Ap [ ], scs_int Ai [ ], scs_float Ax [ ],
    scs_int Lp [ ], scs_int Parent [ ], scs_int Lnz [ ], scs_int Li [ ], scs_float Lx [ ],
    scs_float D [ ], scs_float Y [ ], scs_int Pattern [ ], scs_int Flag [ ],
    scs_int P [ ], scs_int Pinv [ ]) ;

void ldl_lsolve (scs_int n, scs_float X [ ], scs_int Lp [ ], scs_int Li [ ],
    scs_float Lx [ ]) ;

void ldl_dsolve (scs_int n, scs_float X [ ], scs_float D [ ]) ;

void ldl_ltsolve (scs_int n, scs_float X [ ], scs_int Lp [ ], scs_int Li [ ],
    scs_float Lx [ ]) ;

void ldl_perm  (scs_int n, scs_float X [ ], scs_float B [ ], scs_int P [ ]) ;
void ldl_permt (scs_int n, scs_float X [ ], scs_float B [ ], scs_int P [ ]) ;

scs_int ldl_valid_perm (scs_int n, scs_int P [ ], scs_int Flag [ ]) ;
scs_int ldl_valid_matrix ( scs_int n, scs_int Ap [ ], scs_int Ai [ ]) ;

/* ========================================================================== */
/* === LDL version ========================================================== */
/* ========================================================================== */

#define LDL_DATE "May 4, 2016"
#define LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define LDL_MAIN_VERSION 2
#define LDL_SUB_VERSION 2
#define LDL_SUBSUB_VERSION 6
#define LDL_VERSION LDL_VERSION_CODE(LDL_MAIN_VERSION,LDL_SUB_VERSION)

