/* 
 * File:   test_broyden.h
 * Author: Pantelis Sopasakis
 *
 * Created on May 26, 2017, 3:33 PM
 */

#ifndef TEST_BROYDEN_H
#define TEST_BROYDEN_H

#include "unit_test_util.h"
#include "glbopts.h"
#include "directions.h"
#include "scs.h"

#ifdef __cplusplus
extern "C" {
#endif    
    
    bool test_cache_increments(char** str);
    
    bool test_nonfull_cache(char** str);
    
    bool test_full_cache(char** str);


#ifdef __cplusplus
}
#endif

#endif /* TEST_BROYDEN_H */

