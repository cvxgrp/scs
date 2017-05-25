#ifndef DIRECTIONS_H
#define DIRECTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

/** 
 * The cache is full-memory
 */    
#define YS_CACHE_FULL 3
/**
 * The cache has been incremented
 */    
#define YS_CACHE_INCREMENT 2
/** 
 * The head of the cache has been reset to 0
 */
#define YS_CACHE_RESET 1
    
    /**
     * Resets the cache. This methods does not free the memory allocated by the 
     * cache, nor does it overwrite the previously cached data. It simply sets the 
     * parameters `mem_current` and `mem_idx` to 0.
     * @param ys_cache
     * @return status code
     */
    scs_int resetYSCache(SUCache * ys_cache);

    /**
     * Restarted Broyden (as it is reported in the paper)
     * @param w
     * @param iter
     * @return status code
     */
    scs_int computeLSBroyden(Work *work, SUCache *ys_cache, scs_int iter);




#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

