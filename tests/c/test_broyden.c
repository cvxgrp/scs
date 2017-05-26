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
    work->stgs->thetabar = 0.8;
    work->su_cache = scs_calloc(1, sizeof (SUCache));
    work->su_cache->S = scs_malloc((1 + memory) * l_size * sizeof (scs_float));
    work->su_cache->U = scs_malloc((1 + memory) * l_size * sizeof (scs_float));
    work->su_cache->mem = memory;
    work->su_cache->mem_current = 0;
    work->Sk = scs_malloc(l_size * sizeof (scs_float));
    work->Yk = scs_malloc(l_size * sizeof (scs_float));
    work->dir = scs_malloc(l_size * sizeof (scs_float));
    work->R = scs_malloc(l_size * sizeof (scs_float));
    work->stepsize = 1;

    randomize_values(work->Sk, l_size);
    randomize_values(work->Yk, l_size);
    randomize_values(work->dir, l_size);
    randomize_values(work->R, l_size);
    randomize_values(work->su_cache->S, l_size * (1 + memory));
    randomize_values(work->su_cache->U, l_size * (1 + memory));
}

static void destroy_work(Work * work) {
    free(work->stgs);
    work->stgs = NULL;
    free(work->Sk);
    work->Sk = NULL;
    free(work->Yk);
    work->Yk = NULL;
    free(work->R);
    work->R = NULL;
    free(work->su_cache->S);
    work->su_cache->S = NULL;
    free(work->su_cache->U);
    work->su_cache->U = NULL;
    free(work->su_cache);
    work->su_cache = NULL;
    free(work);
    work = NULL;
}

bool test_cache_increments(char** str) {
    Work * work = scs_calloc(1, sizeof (Work));
    scs_int i;
    scs_int method_status;
    prepare_work(work, 3, 10);
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_current, 0, str, "initial mem not 0")
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem, 10, str, "memory not set correctly")

    for (i = 0; i < 100; ++i) {
        method_status = computeLSBroyden(work);
        if (i > 1 && work->su_cache->mem_current == 0) {
            ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_RESET, str, "status not SU_CACHE_RESET")
        } else {
            ASSERT_EQAUL_INT_OR_FAIL(method_status, SU_CACHE_INCREMENT, str, "status not SU_CACHE_INCREMENT")
        }
        ASSERT_TRUE_OR_FAIL(work->su_cache->mem_current <= work->su_cache->mem, 
                str, "mem of cache overflowed")
    }

    resetSUCache(work->su_cache);
    ASSERT_EQAUL_INT_OR_FAIL(work->su_cache->mem_current, 0, str, "active mem after reset not 0")

    destroy_work(work);

    SUCCEED(str)
}

bool test_nonfull_cache(char** str) {
    SUCCEED(str)
}

bool test_full_cache(char** str) {
    SUCCEED(str)
}