#include "directions.h"


scs_int resetYSCache(SUCache * ys_cache) {
    ys_cache->mem_current = 0;
    ys_cache->mem_idx = -1;
    return YS_CACHE_RESET;
}

scs_int computeLSBroyden(Work *work, SUCache *cache, scs_int iter) {
    /* declarations */
    scs_int i;
    
    /* d [work->dir] = -R [work->R] */
    setAsScaledArray(work->dir, work->R, -1.0, work->l);
    
    /* s_tilde_current = y [work->Yk] */
    
    /* update s and d */
    
    
    return 0;
}