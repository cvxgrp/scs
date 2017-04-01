/* 
 * File:   unitTests.c
 * Author: Chung
 *
 * Created on April 1, 2017, 2:10 AM
 */

#include "unitTests.h"

int testOK(char** str) {
    *str = (char*) MESSAGE_OK;
    return TEST_SUCCESS;
}

int testFail(char** str) {
    FAIL_WITH_MESSAGE(str, "Assertion Failed");
}

int testJustAnother(char** str) {
    scs_float X[5];
    scs_float Y[5];
    unsigned int i;

    /* test whether two integers are equal */
    if (!assertEqualsInt(10, 10)) {
        FAIL_WITH_MESSAGE(str, "Integers not equal");
    }
    /* test whether two float are almost equal */
    if (!assertEqualsFloat(10.0, 10.0 + 1e-6, 1e-5)) {
        FAIL_WITH_MESSAGE(str, "Float not equal");
    }

    for (i = 0; i < 5; i++) {
        X[i] = 5.523432 * (i + 1);
        Y[i] = X[i] + 1e-6;
    }

    if (!assertEqualsArray(X, Y, 5, 1e-5)) {
        FAIL_WITH_MESSAGE(str, "The two arrays are not equal");
    }
    
    *str = (char *) MESSAGE_OK;
    return TEST_SUCCESS;
}

int main(int argc, char** argv) {
    int r = TEST_SUCCESS;
    r += test(&testOK, "Dummy passing test");
    r += test(&testFail, "Doomed to fail");
    r += test(&testJustAnother, "Just another test");
    if (r == TEST_SUCCESS) {
        return (EXIT_SUCCESS);
    } else {
        return (EXIT_FAILURE);  
    }
}

