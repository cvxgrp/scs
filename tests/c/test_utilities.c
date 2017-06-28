#include "test_utilities.h"
#include <time.h>

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

    SUCCEED(str);
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
    SUCCEED(str); /* if it reaches this point, it has succeeded */
}

bool testGemm(char** str) {
    const double A[6] = {0.8147,
        0.9058,
        0.1270,
        0.9134,
        0.6324,
        0.0975};
    const double B[6] = {0.2785,
        0.5469,
        0.9575,
        0.9649,
        0.1576,
        0.9706};
    double C[4] = {0.9572,
        0.4854,
        0.8003,
        0.1419};
    const double Cexp[4] = {2.3653,
        1.3934,
        2.3106,
        0.8401};
    const double alpha = 0.5;
    const double beta = 2;
    dgemm_nn(2, 2, 3,
            alpha, A,
            1, 2,
            B,
            1, 3,
            beta, C,
            1, 2);

    ASSERT_EQUAL_ARRAY_OR_FAIL(C, Cexp, 4, 1e-4, str, "gemm failed");
    SUCCEED(str);

}

bool testGemmCP(char** str) {
    const double A[10] = {
        0.334430155748757,
        -0.119893174350795,
        0.804005060428243,
        -0.067975658376914,
        -1.063394117875069,
        0.809765549484799,
        -1.878509454657401,
        -0.259736781357468,
        0.470502094834760,
        0.320050300435137
    };
    const double B[15] = {
        0.242754814623263,
        -0.103589012251697,
        -0.454961543295210,
        -0.413269912824790,
        1.497631598138995,
        -0.128084953472689,
        2.266925189882620,
        -0.254500560127930,
        -0.711248533058385,
        -0.369507289387400,
        -1.999207676978967,
        -0.555866380148284,
        0.587186741401126,
        1.004782851967581,
        -0.585280006989040
    };
    double C[6] = {
        -0.774073198694521,
        -0.960044604299499,
        -2.022434124346632,
        -0.323079592516537,
        1.037274430734777,
        0.420892205865074
    };
    const double alpha = -0.286281752586377;
    const double beta = 3.194915595797473;
    const double Cexp[6] = {
        -3.034975746827981,
        -3.123425247115062,
        -7.381229796362662,
        -0.952525926453145,
        4.431303305975694,
        1.257495836652682
    };

    matrixMultiplicationColumnPacked(2, 3, 5, alpha, A, beta, B, C);

    ASSERT_EQUAL_ARRAY_OR_FAIL(C, Cexp, 4, 1e-14, str, "gemm failed");

    SUCCEED(str);
}

bool testGemmTransCP(char** str) {
    const double A[12] = {
        0.698299663682011,
        -1.627423017907931,
        -1.372695305499414,
        -1.100828998920425,
        1.619000819707846,
        -0.600157916750174,
        -0.540089717235530,
        1.484871682894813,
        1.809840858337382,
        0.919301984685824,
        -0.212130772097334,
        -0.095040503915385
    };
    const double B[8] = {
        0.701256481812284,
        0.876974554050047,
        -2.190732553342963,
        0.687223989397896,
        0.905368244420720,
        2.186309802484150,
        -0.496517337448137,
        0.288763931098904
    };
    double C[6] = {
        -1.608876042935446,
        -0.040192422065262,
        1.723531705742089,
        0.445855130092155,
        -0.628575736932150,
        -0.462395267263025
    };
    const double alpha = -0.023912990352431;
    const double beta = 0.916952300228893;
    const double Cexp[6] = {
        -1.506664428673252,
        -0.104113242719988,
        1.521217097262638,
        0.470096441685509,
        -0.596714407327636,
        -0.513102186175089
    };
    matrixMultiplicationTransColumnPacked(3, 2, 4, alpha, A, beta, B, C);
    ASSERT_EQUAL_ARRAY_OR_FAIL(C, Cexp, 4, 1e-5, str, "gemm failed");
    SUCCEED(str);
}

bool testUnrolledDot(char** str) {
    const scs_float x[4] = {5., 2., 3., 4.};
    const scs_float y[4] = {0.1, 0.2, 0.5, -0.9};
    const scs_float a[5] = {5.0, 6.0, 7.0, 8.0, 9.0};
    const scs_float b[5] = {11.0, 13.0, 15.0, 17.0, 19.0};
    const scs_float c[6] = {-1.2, 3.4, 5.6, 7.8, 9.1, 11.12};
    const scs_float d[6] = {2.3, 4.5, 6.7, 8.9, 10.11, 12.13};
    const scs_float p[11] = {0.554259561551264, -2.159803537045049, -0.372961114202917,
        -2.592374718571431, 0.549237254991523, 0.535999299756211, -0.188665249567086,
        0.760181191234408, 1.033831426574950, -0.196086216446758, 0.474596953309341};
    const scs_float q[11] = {-1.2360603341291734, -0.7844101929878990, 0.7651738841267504,
        -0.8603726595864287, 1.0311516752375671, 0.5283233936417885, -1.2091816279319294,
        -0.1549870848290139, -0.9608099563757425, 0.7599149377251431, 0.0224256807626837};
    /*
        const unsigned n = 5e6;
        scs_float * big_x = malloc(n * sizeof (scs_float));
        scs_float * big_y = malloc(n * sizeof (scs_float));    
        scs_int i;
     */
    
    scs_float ip;
    /*
        clock_t t;
     */

    ip = innerProd(x, y, 4);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, -1.2, 1e-6, str, "inn prod 1");
    ip = innerProd(a, b, 5);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 545.0, 1e-6, str, "inn prod 2");
    ip = innerProd(c, d, 6);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 346.3666, 1e-6, str, "inn prod 3");
    ip = innerProd(p, q, 11);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 2.78226083574042, 1e-10, str, "inn prod 4");
    ip = innerProd(p, q, 7);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 4.03176031557307, 1e-10, str, "inn prod 5");
    ip = innerProd(p, q, 8);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 3.91394204880181, 1e-10, str, "inn prod 6");
    ip = innerProd(p, q, 1);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, -0.685098258845344, 1e-10, str, "inn prod 7");
    ip = innerProd(p, q, 2);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 1.00907365046411, 1e-10, str, "inn prod 8");
    ip = innerProd(p, q, 0);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 0.0f, 1e-10, str, "inn prod 9");
    ip = innerProd(p, q, -1);
    ASSERT_EQAUL_FLOAT_OR_FAIL(ip, 0.0f, 1e-10, str, "inn prod 10");

    /*
        t = clock();
        for (i = 0; i < 500; ++i) {
            ip = innerProd(big_x, big_y, n);
        }
        t = clock() - t;
        double elapsed_time = t / (double) CLOCKS_PER_SEC;
        printf("Elapsed time: %g.\n", elapsed_time);
        scs_free(big_x);
        scs_free(big_y);
     */
    SUCCEED(str);
}

bool testLinAlg(char** str) {
#define n 97
    scs_float x[n];
    scs_float y[n];
    unsigned int i = 0;
    for (i = 0; i < n; ++i) {
        x[i] = SQRTF(0.1*i+1);
        y[i] = sin(i/20);
    }
    subtractArray(x, y, n);
    for (i = 0; i < n; ++i) {
        ASSERT_EQAUL_FLOAT_OR_FAIL(x[i], SQRTF(0.1*i+1) - sin(i/20), 1e-10, str, "wrong");
    }
    SUCCEED(str);
}