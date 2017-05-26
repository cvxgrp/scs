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
     * iteration (current FPR, values of S and Y).
     * 
     * @return status code of the method.
     */
    scs_int computeLSBroyden(Work *work);


#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

