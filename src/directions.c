#include "directions.h"

scs_int pushToYSCache(Work *work, YSCache * ys_cache) {

    scs_int cache_status = YS_CACHE_INCREMENT;

    /* Declarations */
    scs_float *Y;
    scs_float *S;
    scs_float *YS;
    scs_float *yk;
    scs_float *sk;

    if (ys_cache->mem_current == ys_cache->mem) {
        cache_status = YS_CACHE_FULL;
    }

    /* increment the memory (when initialised, mem_idx = -1) */
    ys_cache -> mem_idx++;
    if (ys_cache -> mem_idx >= ys_cache->mem) {
        cache_status = YS_CACHE_RESET;
        ys_cache -> mem_idx = 0;
    }

    /* update the current memory */
    ys_cache->mem_current++; /* increase current memory */
    if (ys_cache->mem_current > ys_cache->mem) {
        ys_cache->mem_current = ys_cache->mem; /* buffer is now full */
    }

    /* push Yk into the right position in the cache */
    Y = ys_cache->Y;
    yk = work->Yk;
    Y += (ys_cache -> mem_idx) * work->l;
    memcpy(Y, yk, work->l * sizeof (scs_float));

    /* push Sk into the right position in the cache */
    S = ys_cache->S;
    sk = work->Sk;
    S += (ys_cache -> mem_idx) * work->l;
    memcpy(S, sk, work->l * sizeof (scs_float));

    /* Compute the inner product and store in the right position in the cache */
    YS = ys_cache->YS;
    YS[ys_cache -> mem_idx] = innerProd(work->Yk, work->Sk, work->l);

    return cache_status;
}

scs_int resetYSCache(YSCache * ys_cache) {
    ys_cache->mem_current = 0;
    ys_cache->mem_idx = -1;
    return YS_CACHE_RESET;
}

void computeLBroyden(Work *w, scs_int iter) {

}

void computeLMBroyden(Work *w, scs_int iter) {


}