#include "test_broyden.h"

static void randomize_values(scs_float* x, scs_int l) {
    scs_int i;
    for (i = 0; i < l; ++i) {
        x[i] = (i + 0.5) * 0.05;
    }
}

static void prepare_work(Work * work, scs_int l_size, scs_int memory) {
    work->l = l_size;
    work->stgs = scs_calloc(1, sizeof (Settings));
    work->stgs->thetabar = 0.2;
    work->su_cache = scs_calloc(1, sizeof (SUCache));
    work->su_cache->S = scs_calloc((1 + memory) * l_size, sizeof (scs_float));
    work->su_cache->U = scs_calloc((1 + memory) * l_size, sizeof (scs_float));
    work->su_cache->mem = memory;
    work->su_cache->mem_cursor = 0;
    work->Sk = scs_calloc(l_size, sizeof (scs_float)); /* malloc would be just fine anyway... */
    work->Yk = scs_calloc(l_size, sizeof (scs_float));
    work->dir = scs_calloc(l_size, sizeof (scs_float));
    work->R = scs_calloc(l_size, sizeof (scs_float));
    work->stepsize = 1;

    randomize_values(work->Sk, l_size);
    randomize_values(work->Yk, l_size);
    randomize_values(work->dir, l_size);
    randomize_values(work->R, l_size);
}

static void destroy_work(Work * work) {
    if (!work) {
        return;
    }
    if (work->stgs) {
        free(work->stgs);
        work->stgs = NULL;
    }
    if (work->Sk) {
        free(work->Sk);
        work->Sk = NULL;
    }
    if (work->Yk) {
        free(work->Yk);
        work->Yk = NULL;
    }
    if (work->R) {
        free(work->R);
        work->R = NULL;
    }
    if (work->dir) {
        free(work->dir);
        work->dir = NULL;
    }
    if (work->su_cache) {
        if (work->su_cache->S) {
            free(work->su_cache->S);
            work->su_cache->S = NULL;
        }
        if (work->su_cache->U) {
            free(work->su_cache->U);
            work->su_cache->U = NULL;
        }
        free(work->su_cache);
        work->su_cache = NULL;
    }
    free(work);
    work = NULL;
}

bool test_cache_increments(char** str) {
    Work * work = scs_calloc(1, sizeof (Work));
    scs_int i;
    const scs_int l = 4;
    const scs_int mem = 10;
    const scs_int runs = 5000;
    scs_int method_status;
    prepare_work(work, l, mem);
    resetSUCache(work->su_cache);
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_cursor, 0, str, "active mem after reset not 0");

    ASSERT_EQAUL_INT_OR_FAIL(work->l, l, str, "size not correct")
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem, mem, str, "wrong memory")
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_cursor, 0, str, "initial mem not 0")
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem, mem, str, "memory not set correctly")

    for (i = 0; i < runs; ++i) {
        method_status = computeLSBroyden(work);
        if (i > 1 && work->su_cache->mem_cursor == 0) {
            ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_RESET, str, "status not SU_CACHE_RESET")
        } else {
            ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_INCREMENT, str, "status not SU_CACHE_INCREMENT")
        }
        ASSERT_TRUE_OR_FAIL(work->su_cache->mem_cursor <= work->su_cache->mem,
                str, "mem of cache overflowed")
    }

    resetSUCache(work->su_cache);
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_cursor, 0, str, "active mem after reset not 0");

    destroy_work(work);

    SUCCEED(str);
}

bool test_broyden_direction_empty_memory(char** str) {
    Work * work = scs_calloc(1, sizeof (Work));
    const scs_int l = 4;
    const scs_int mem = 10;
    scs_int method_status;

    prepare_work(work, l, mem);

    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem, mem, str, "memory not set correctly");
    work->stgs->thetabar = 0.4;

    /* Yk = [1.5; 3.2; -1.8; 0.7]; */
    work->Yk[0] = 1.5;
    work->Yk[1] = 3.2;
    work->Yk[2] = -1.8;
    work->Yk[3] = 0.7;

    /* Sk = [0.5; 4.1; -3.0; 6.2]; */
    work->Sk[0] = 0.5;
    work->Sk[1] = 4.1;
    work->Sk[2] = -3.0;
    work->Sk[3] = 6.2;

    /* R = [3.3; 5.1, 0.7, -5.2];*/
    work->R[0] = 3.3;
    work->R[1] = 5.1;
    work->R[2] = 0.7;
    work->R[3] = -5.2;

    work->stepsize = 0.9;

    method_status = computeLSBroyden(work);
    ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_INCREMENT, str, "memory not incremented");

    scs_float u_expected[4] = {-0.0366837857666911, 0.0330154071900220, -0.0440205429200293, 0.2017608217168012};
    scs_float d_expected[4] = {-3.73213499633162, -4.71107850330154, -1.21856199559795, 7.57674247982392};
    ASSERT_EQUAL_ARRAY_OR_FAIL(work->su_cache->U, u_expected, l, 1e-10, str, "u not correct");
    ASSERT_EQUAL_ARRAY_OR_FAIL(work->dir, d_expected, l, 1e-10, str, "direction not correct");
    ASSERT_EQUAL_ARRAY_OR_FAIL(work->su_cache->S, work->Sk, l, 1e-10, str, "sk not added to the cache");


    destroy_work(work);

    SUCCEED(str);
}

bool test_cache_s(char** str) {
    Work * work = scs_calloc(1, sizeof (Work));
    scs_int i;
    scs_int j;
    const scs_int l = 4;
    const scs_int mem = 5;
    scs_int method_status;

    prepare_work(work, l, mem);


    for (i = 0; i < 20; ++i) {

        work->Sk[0] = 1.0 + 4.5 / i;
        work->Sk[1] = 50.1 + 16.2 / i;
        work->Sk[2] = 10.0 / i;
        work->Sk[3] = 1.0 / i;


        work->R[0] *= 0.01;
        work->R[1] *= 0.01;
        work->R[2] *= 0.01;
        work->R[3] *= 0.01;

        work->Yk[0] *= 0.1;
        work->Yk[1] *= 0.1;
        work->Yk[2] *= 0.1;
        work->Yk[3] *= 0.1;

        method_status = computeLSBroyden(work);


        if ((i + 1) % (mem) == 0) {
            ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_cursor, 0, str, "current mem not zero");
            ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_RESET, str, "not reset");
        } else {
            ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_INCREMENT, str, "not reset");
            ASSERT_TRUE_OR_FAIL(work->su_cache->mem_cursor > 0, str, "memory cursor is at zero");
            ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_cursor, (i + 1) % (mem), str, "cursor at wrong position")
        }
    }

    destroy_work(work);

    SUCCEED(str);
}