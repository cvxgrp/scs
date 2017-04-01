/* 
 * File:   unitTests.h
 * Author: chung
 *
 * Created on April 1, 2017, 2:14 AM
 */

#ifndef UNITTESTS_H
#define UNITTESTS_H

#include <stdio.h>
#include <stdlib.h>
#include "glbopts.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TEST_SUCCESS 0 /**< test is successful */
#define TEST_FAILURE 1 /**< test fails */
#define MESSAGE_OK "OK" /**< a message returned when a test is successful */
#define MESSAGE_MAX_LENGTH 1024 /**< maximum message length */

    /**
     * Fails with a given message.
     * FAIL_WITH_MESSAGE(str, message)
     */
#define FAIL_WITH_MESSAGE(str, message)\
        *str = (char*)message;\
        return TEST_FAILURE;

    /**
     * Function template defining a unit test:
     * 
     *  int myTestFunction(char**);
     * 
     * This type is a pointer to such a function which takes as an input argument 
     * a pointer to a string (char**) and returns a status code (either TEST_SUCCESS
     * or TEST_FAILURE).
     */
    typedef int (*unitTest_t)(char**);

    /**
     * Dummy successful test.
     * This serves only as an example.
     * @param msg message ("OK")
     * @return returns #TEST_SUCCESS.
     */
    int testOK(char** msg);

    /**
     * Tester function.
     * @param ut Unit Test function handle
     * @param name Name of the test
     * @return TEST_SUCCESS if the test succeeds and TEST_FAILURE if it fails.
     */
    int test(const unitTest_t ut, const char* name);

    /**
     * Assert that two integers are equal.
     * @param a
     * @param b
     * @return 
     */
    int assertEqualsInt(const scs_int a, const scs_int b);

    /**
     * Assert that two floats are equal up to a given tolerance.
     * @param a
     * @param b
     * @param tol tolerance
     * @return 
     */
    int assertEqualsFloat(const scs_float a, const scs_float b, const scs_float tol);

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
    int assertEqualsArray(
            const scs_float * a,
            const scs_float * b,
            scs_int n,
            const scs_float tol);

    /* ---------- SOME IMPLEMENTATIONS --------------- */

    int assertEqualsInt(const scs_int a, const scs_int b) {
        return (a == b);
    }

    int assertEqualsFloat(const scs_float a, const scs_float b, const scs_float tol) {
        return ( fabs(a - b) < tol);
    }

    int assertEqualsArray(
            const scs_float * a,
            const scs_float * b,
            scs_int n,
            const scs_float tol) {
        unsigned int i;
        for (i = 0; i < n; i++) {
            if (!assertEqualsFloat(a[i], b[i], tol)) {
                return 0;
            }
        }
        return 1;
    }

    int test(const unitTest_t ut, const char* name) {
        char * message = malloc(MESSAGE_MAX_LENGTH * sizeof (char));
        int result = ut(&message);
        if (result == TEST_SUCCESS) {
            printf("[ PASS ]");
        } else {
            printf("< FAIL >");
        }
        printf(" (%s) -- %s\n", name, message);
        return result;
    }


#ifdef __cplusplus
}
#endif

#endif /* UNITTESTS_H */

