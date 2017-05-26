#ifndef DIRECTIONS_H
#define DIRECTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

/** 
 * The cache is full-memory
 */    
#define SU_CACHE_FULL 3
/**
 * The cache has been incremented
 */    
#define SU_CACHE_INCREMENT 2
/** 
 * The head of the cache has been reset to 0
 */
#define SU_CACHE_RESET 1
    
    /**
     * Resets the cache. This methods does not free the memory allocated by the 
     * cache, nor does it overwrite the previously cached data. It simply sets the 
     * parameters `mem_current` and `mem_idx` to 0.
     * @param cache
     * @return status code
     */
    scs_int resetSUCache(SUCache * cache);

    /**
     * Restarted Broyden (as it is reported in the paper).
     * 
     * @param work Work structure with all available information about the current
     * iteration (current FPR, values of S and Y). Work must provide the following 
     * information
     * 
     *   - <code>work->stgs->thetabar</code>
     *   - <code>work->su_cache</code> with
     *      - <code>work->su_cache->S</code> 
     *      - <code>work->su_cache->U</code> 
     *      - <code>work->su_cache->current_mem</code> 
     *      - <code>work->su_cache->mem</code> (nonzero memory)
     *   - <code>work->Yk</code>
     *   - <code>work->Sk</code>
     *   - <code>work->R</code>
     *   - <code>work->stepsize</code>
     *   - <code>work->l</code>
     * 
     * @return status code of the method.
     */
    scs_int computeLSBroyden(Work *work);


#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

