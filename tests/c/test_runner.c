
#include "test_runner.h"
#include "test_superscs.h"


int main(int argc, char** argv) {
    int r;
    number_of_assertions = 0;    
    
    printf("\n***** Test Results *****\n\n");    
    
    r = TEST_SUCCESS;
    
    /* Test functions: */   
    r += test(&test_dummy_method, "Dummy passing test");
    r += test(&test_superscs_011_progress, "Test SuperSCS (0,1,1) thoroughly");
    r += test(&testMillisToTime, "Milliseconds to time");
    r += test(&testUnrolledDot, "Unrolled dot");
    r += test(&testLinAlg, "Unrolled subtraction");
    r += test(&testProjLinSysv2, "Test projLinSysv2");
    r += test(&testScaleArray, "Test scaleArray");
    r += test(&testGemm, "Test GEMM");
    r += test(&testGemmCP, "Test GEMM/CP");
    r += test(&testGemmTransCP, "Test GEMM/Tran/CP");
    r += test(&test_cache_increments, "Test Restarted Broyden cache");
    r += test(&test_broyden_direction_empty_memory, "Test Restarted Broyden direction");     
    r += test(&test_full_broyden, "Test full Broyden");     
    r += test(&test_cache_s, "Test Broyden S-Cache");
    r += test(&test_broyden, "Test Broyden dir correctness");
    r += test(&test_superscs_solve, "Test SuperSCS");
    r += test(&test_superscs_000, "Test SuperSCS (0,0,0)");
    r += test(&test_superscs_001_fpr, "Test SuperSCS (0,0,1) with FPR"); 
    r += test(&test_superscs_001_rbroyden, "Test SuperSCS (0,0,1) with R-Broyden");
    r += test(&test_superscs_100_rbroyden, "Test SuperSCS (1,0,0) with R-Broyden"); 
    r += test(&test_residuals, "Test residuals"); 
    r += test(&test_rho_x, "Test rho_x");
    r += test(&test_validation, "Test validation");
    r += test(&test_no_normalization, "Test SuperSCS unnormalized");     
    r += test(&test_warm_start, "Test SuperSCS warm_start");
    r += test(&test_scale, "Test SuperSCS scalings");
    printf("\nTotal assertions: %d\n", number_of_assertions);
    if (r == TEST_SUCCESS) {
        printf("\n~ All tests passed\n\n");
        return (EXIT_SUCCESS);
    } else {
        printf("\n~ %d Tests Failed\n\n", r);
        return (EXIT_FAILURE);
    }
    
}