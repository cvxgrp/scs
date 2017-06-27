/* 
 * File:   test_utilities.h
 * Author: Pantelis Sopasakis
 *
 * Created on May 26, 2017, 3:25 PM
 */

#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include "unit_test_util.h"

#ifdef __cplusplus
extern "C" {
#endif

    bool testProjLinSysv2(char** str);

    bool testScaleArray(char** str);
    
    bool testGemm(char** str);
    
    bool testGemmCP(char** str);
    
    bool testGemmTransCP(char** str);
    
    bool testUnrolledDot(char** str);


#ifdef __cplusplus
}
#endif

#endif /* TEST_UTILITIES_H */

