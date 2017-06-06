#ifndef DIRECTIONS_H
#define DIRECTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

/**
 * The cache has been incremented.
 */    
#define SU_CACHE_INCREMENT 2
/** 
 * The cursor of the cache has been reset to \c 0.
 */
#define SU_CACHE_RESET 1
    
    /**
     * Resets the cache. This methods does not free the memory allocated by the 
     * cache, nor does it overwrite the previously cached data. It simply sets the 
     * parameters `mem_current` and `mem_idx` to 0.
     * 
     * @param cache the cache to be reset 
     * 
     * @return status code (returns #SU_CACHE_RESET)
     */
    scs_int resetSUCache(SUCache * cache);

    /**
     * Restarted Broyden (as it is reported in the paper).
     * 
     * @param work Work structure with all available information about the current
     * iteration (current FPR, values of \f$s_k\f$, \f$y_k\f$ etc).
     * 
     * @return status code of the method.
     * 
     * @see \ref sec-restarted-broyden "Restarted Broyden Algorithm"
     */
    scs_int computeLSBroyden(Work *work);
    
    /**
     * Computes a direction according to the value of 
     * <code>work->stgs->direction</code>.
     * 
     * @param work
     * @return 
     * 
     * @see direction_type
     */
    scs_int computeDirection(Work *work);


#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

