#include "test_utilities.h"

bool testProjLinSysv2(char** str) {

    scs_int n = 5, m = 10, l = n + m + 1, i;
    scs_float *u_t, *u, *h, *g;
    scs_float gTh = 2.2;
    const scs_float rho_x = 1.0;
    scs_float *expected_result;
    bool test_pass = false;
    float tolerance = 1e-6;


    /* memory allocation */
    expected_result = malloc(l * sizeof (scs_float));
    u = malloc(l * sizeof (scs_float));
    u_t = malloc(l * sizeof (scs_float));
    h = malloc((l - 1) * sizeof (scs_float));
    g = malloc((l - 1) * sizeof (scs_float));


    for (i = 0; i < l; ++i) {
        u_t[i] = 0.5 * (i + 1);
    }

    for (i = 0; i < l - 1; ++i) {
        h[i] = 0.2 * (i + 1);
        g[i] = 0.8 * (i + 1);
    }
    /*  memcpy(u_t, u, l * sizeof (scs_float)); */

    scaleArray(u_t, rho_x, n);

    addScaledArray(u_t, h, l - 1, -u_t[l - 1]);
    addScaledArray(u_t, h, l - 1,
            -innerProd(u_t, g, l - 1) / (gTh + 1));
    scaleArray(&(u_t[n]), -1, m);
    /*   status = solveLinSys(A, stgs, w->p, w->u_t, w->u, iter); */
    u_t[l - 1] += innerProd(u_t, h, l - 1);

    expected_result[0] = 67.10;
    expected_result[1] = 134.20;
    expected_result[2] = 201.30;
    expected_result[3] = 268.40;
    expected_result[4] = 335.50;
    expected_result[5] = -402.60;
    expected_result[6] = -469.70;
    expected_result[7] = -536.80;
    expected_result[8] = -603.90;
    expected_result[9] = -671.00;
    expected_result[10] = -738.10;
    expected_result[11] = -805.20;
    expected_result[12] = -872.30;
    expected_result[13] = -939.40;
    expected_result[14] = -1006.50;
    expected_result[15] = -15156.60;


    test_pass = assertEqualsArray(u_t, expected_result, l, tolerance);

    /* free memory */
    free(u);
    u = NULL;
    free(u_t);
    u_t = NULL;
    free(h);
    h = NULL;
    free(g);
    g = NULL;
    free(expected_result);
    expected_result = NULL;
    if (!test_pass) {
        FAIL_WITH_MESSAGE(str, "testProjLinSysv2 failed");
    }

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
    a = NULL;
    free(expected_result);
    expected_result = NULL;
    if (!test_pass) {
        FAIL_WITH_MESSAGE(str, "scaleArray failed");
    }
    SUCCEED /* if it reaches this point, it has succeeded */
}

