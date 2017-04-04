/* 
 * File:   unitTests.c
 * Author: Pantelis Sopasakis
 *
 * Created on April 1, 2017, 2:10 AM
 */

#include "unitTests.h"

bool testOK(char** str) {
    SUCCEED
}

bool testLinsysV2(char** str) {
    SUCCEED
}

bool testScaleArray(char** str) {
    const scs_int N = 10;
    float tolerance = 1e-6;
    unsigned int i;
    bool test_pass = false;
    scs_float * a = malloc(N * sizeof (scs_float));
    scs_float * expected_result = malloc(N * sizeof (scs_float));
    const scs_float b = 3.23412;

    for (i = 0; i < N; ++i) {
        a[i] = 0.5 * (i + 1);
        expected_result[i] = b * a[i];
    }
    scaleArray(a, b, N);
    test_pass = assertEqualsArray(a, expected_result, N, tolerance);

    /* free memory */
    free(a);
    free(expected_result);
    if (!test_pass) FAIL_WITH_MESSAGE(str, "scaleArray failed");
    SUCCEED /* if it reaches this point, it has succeeded */
}

int main(int argc, char** argv) {
    int r = TEST_SUCCESS;
    /* Test functions: */
    r += test(&testOK, "Dummy passing test");
    r += test(&testScaleArray, "Scale array");
    r += test(&testLinsysV2, "test solve lin sys v2");
    if (r == TEST_SUCCESS) {
        return (EXIT_SUCCESS);
    } else {
        return (EXIT_FAILURE);
    }
}

