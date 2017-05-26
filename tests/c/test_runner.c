
#include "test_runner.h"


int main(int argc, char** argv) {
    int r;
    
    printf("Test Results:\n\n");    
    
    r = TEST_SUCCESS;
    /* Test functions: */
    r += test(&test_dummy_method, "Dummy passing test");
    r += test(&testProjLinSysv2, "Test projLinSysv2");
    r += test(&testScaleArray, "Test scaleArray");
    r += test(&test_cache_increments, "Test Broyden cache");
    if (r == TEST_SUCCESS) {
        printf("\nAll tests passed\n");
        return (EXIT_SUCCESS);
    } else {
        printf("\n%d Tests Failed\n", r);
        return (EXIT_FAILURE);
    }
}