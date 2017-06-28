/* 
 * File:   test_superscs.h
 * Author: chung
 *
 * Created on May 31, 2017, 5:41 PM
 */

#ifndef TEST_SUPERSCS_H
#define TEST_SUPERSCS_H

#include "unit_test_util.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    bool test_superscs_solve(char** str);
    
    bool test_superscs_000(char** str);
    
    bool test_superscs_001_fpr(char** str);
    
    bool test_superscs_001_rbroyden(char** str);
    
    bool test_superscs_100_rbroyden(char** str);
    
    bool test_residuals(char** str);

#ifdef __cplusplus
}
#endif

#endif /* TEST_SUPERSCS_H */

