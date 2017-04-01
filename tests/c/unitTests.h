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

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define MESSAGE_OK "OK"
#define MESSAGE_MAX_LENGTH 255

    /**
     * Fails with a given message
     */
    #define FAIL_WITH_MESSAGE(str, message)\
        *str = (char*)message;\
        return TEST_FAILURE;

    /**
     * Function template defining a unit test:
     * 
     *  int myTestFunction(char**);
     * 
     * This type is a pointer to such a function.
     */
    typedef int (*unitTest_t)(char**);

    /**
     * Dummy successful test.
     * @param msg message ("OK")
     * @return returns #TEST_SUCCESS.
     */
    int testOK(char** msg);

    /**
     * Tester function.
     * @param ut Unit Test function handle
     * @param name Name of the test
     */
    void test(const unitTest_t ut, const char* name);

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
     * @param tol
     * @return 
     */
    int assertEqualsFloat(const scs_float a, const scs_float b, const scs_float tol);

    
    
    
    
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

    void test(const unitTest_t ut, const char* name) {
        char * message = malloc(MESSAGE_MAX_LENGTH * sizeof (char));
        int result = ut(&message);
        if (result == TEST_SUCCESS) {
            printf("[ PASS ]");
        } else {
            printf("< FAIL >");
        }
        printf(" (%s) -- %s\n", name, message);
    }


#ifdef __cplusplus
}
#endif

#endif /* UNITTESTS_H */

