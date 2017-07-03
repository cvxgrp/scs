/* 
 * File:   unitTests.h
 * Author: Pantelis Sopasakis
 *
 * Created on April 1, 2017, 2:14 AM
 */

#ifndef UNITTESTS_H
#define UNITTESTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "glbopts.h"
#include <math.h>
#include <stdbool.h>
#include "linAlg.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef bool
#undef bool
#endif

#ifdef true
#undef true
#endif

#ifdef false
#undef false
#endif
    
    typedef int bool;
#define true 1
#define false 0


    int number_of_assertions;

#define TEST_SUCCESS 0 /**< test is successful */
#define TEST_FAILURE 1 /**< test fails */
#define MESSAGE_OK "OK" /**< a message returned when a test is successful */
#define TEST_PASS_FLAG "\x1B[92m[PASS]\x1B[39m " /**< flag for PASS */
#define TEST_FAIL_FLAG "\x1B[31m<FAIL>\x1B[39m " /**< flag for FAIL */   

    /**
     * Fails with a given message.
     */
#define FAIL_WITH_MESSAGE(str, message)\
        *str = (char*)(message);\
        return TEST_FAILURE


#define ASSERT_TRUE_OR_FAIL(p, str, message)\
    number_of_assertions++;\
    if (!(p)) { \
        FAIL_WITH_MESSAGE((str), (message));\
    }

    /**
     * Check whether two integers are equal, or fail with a given message.
     */
#define ASSERT_EQUAL_INT_OR_FAIL(val, expected, str, message)\
        number_of_assertions++;\
        if (!assertEqualsInt((val),(expected))) { \
          {\
           char buff[500];\
           char error_msg[100];\
           sprintf(error_msg, "\n\tExpected: %d, Actual %d", expected, val);\
           strcpy(buff, message);\
           strcat(buff, error_msg);\
           FAIL_WITH_MESSAGE((str), (buff)); \
          }\
        }

    /**
     * Check whether two integers are equal, or fail with a given message.
     */
#define ASSERT_EQUAL_FLOAT_OR_FAIL(val, expected, tol, str, message)\
        number_of_assertions++;\
        if (!assertEqualsFloat((val), (expected), (tol))) {\
           char buff[500];\
           char error_msg[100];\
           sprintf(error_msg, "\n\tExpected: %g, Actual %g (tol=%g)", expected, val, tol);\
           strcpy(buff, message);\
           strcat(buff, error_msg);\
           FAIL_WITH_MESSAGE((str), (buff)); \
        }

    /**
     * Check whether two arrays are equal, or fail with a given message.
     */
#define ASSERT_EQUAL_ARRAY_OR_FAIL(val,expected,len,tol,str,message)\
    number_of_assertions++;\
    if (!assertEqualsArray((val),(expected),(len),(tol))){\
      FAIL_WITH_MESSAGE((str), (message));\
    }

    /**
     * Succeed
     */
#define SUCCEED(str)\
        *str = (char*) MESSAGE_OK;\
        return TEST_SUCCESS

    /**
     * Function template defining a unit test:
     * 
     *  int myTestFunction(char**);
     * 
     * This type is a pointer to such a function which takes as an input argument 
     * a pointer to a string (char**) and returns a status code (either TEST_SUCCESS
     * or TEST_FAILURE).
     */
    typedef bool (*unitTest_t)(char**);


    /**
     * Tester function.
     * @param ut Unit Test function handle
     * @param name Name of the test
     * @return TEST_SUCCESS if the test succeeds and TEST_FAILURE if it fails.
     */
    bool test(const unitTest_t ut, const char* name);

    /**
     * Assert that two integers are equal.
     * @param a
     * @param b
     * @return 
     */
    bool assertEqualsInt(const scs_int a, const scs_int b);

    /**
     * Assert that two floats are equal up to a given tolerance.
     * @param a
     * @param b
     * @param tol tolerance
     * @return 
     */
    bool assertEqualsFloat(const scs_float a, const scs_float b, const scs_float tol);

    /**
     * Checks whether two arrays of float are equal, element-wise, up to a certain
     * tolerance.
     * 
     * @param a first array
     * @param b second array
     * @param n length of array
     * @param tol tolerance
     * @return 
     */
    bool assertEqualsArray(
            const scs_float * a,
            const scs_float * b,
            scs_int n,
            const scs_float tol);


#ifdef __cplusplus
}
#endif

#endif /* UNITTESTS_H */

