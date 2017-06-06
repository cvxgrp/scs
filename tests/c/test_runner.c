
#include "test_runner.h"
#include "test_superscs.h"


int main(int argc, char** argv) {
    int r;
    
    printf("\n***** Test Results *****\n\n");    
    
    r = TEST_SUCCESS;
    
    /* Test functions: */
    
    r += test(&test_dummy_method, "Dummy passing test");
    r += test(&testProjLinSysv2, "Test projLinSysv2");
    r += test(&testScaleArray, "Test scaleArray");
    r += test(&test_cache_increments, "Test Broyden cache");
    r += test(&test_broyden_direction_empty_memory, "Test Broyden direction");     
    r += test(&test_cache_s, "Test Broyden S-Cache");
    r += test(&test_broyden, "Test Broyden dir correctness");
    r += test(&test_superscs, "Test SuperSCS");
    
    if (r == TEST_SUCCESS) {
        printf("\n~ All tests passed\n\n");
        return (EXIT_SUCCESS);
    } else {
        printf("\n~ %d Tests Failed\n\n", r);
        return (EXIT_FAILURE);
    }
}