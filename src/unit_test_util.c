/* 
 * File:   unitTests.c
 * Author: Pantelis Sopasakis
 *
 * Created on April 1, 2017, 2:10 AM
 */

#include "unit_test_util.h"
#include "directions.h"

bool assertEqualsInt(const scs_int a, const scs_int b) {
    return (a == b);
}

bool assertEqualsFloat(const scs_float a, const scs_float b, const scs_float tol) {
    return ( fabs(a - b) < tol);
}

bool assertEqualsArray(
        const scs_float * a,
        const scs_float * b,
        scs_int n,
        const scs_float tol) {
    return calcNormInfDiff(a, b, n) < tol;
}

bool test(const unitTest_t ut, const char* name) {
    char * message = NULL;
    int result = ut(&message);
    if (result == TEST_SUCCESS) {
        printf(TEST_PASS_FLAG);
    } else {
        printf(TEST_FAIL_FLAG); /* LCOV_EXCL_LINE */
    }
    printf(" (%s) -- %s\n", name, message);
    return result;
}



