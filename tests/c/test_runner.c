
#include "test_runner.h"


int main(int argc, char** argv) {
    int r = TEST_SUCCESS;
    /* Test functions: */
    r += test(&test_dummy_method, "Dummy passing test");
    r += test(&testProjLinSysv2, "Test projLinSysv2");
    r += test(&testScaleArray, "Test scaleArray");
    if (r == TEST_SUCCESS) {
        return (EXIT_SUCCESS);
    } else {
        return (EXIT_FAILURE);
    }
}