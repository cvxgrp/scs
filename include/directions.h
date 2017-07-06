#ifndef DIRECTIONS_H
#define DIRECTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"

    /**
     * The cache has been incremented.
     */
#define SU_CACHE_INCREMENT 101
    /** 
     * The cursor of the cache has been reset to \c 0.
     */
#define SU_CACHE_RESET 100
    /**
     * The direction could not be computed due to an error.
     */
#define DIRECTION_ERROR -1
    /**
     * The direction was computed successfully.
     * 
     * All nonnegative status codes denote success.
     */
#define DIRECTION_SUCCESS 0

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
     * Full Broyden method.
     * 
     * @param work
     * @param i
     * 
     * @return status code of the method.
     * 
     * @see \ref sec-full-broyden "Full Broyden Algorithm"
     * 
     * \warning Not implemented yet
     */
    scs_int computeFullBroyden(Work *work, scs_int i);

    /**
     * Frees memory allocated for the full Broyden method.
     */
    void freeFullBroyden(void);

    /**
     * Computes a direction according to the value of 
     * <code>work->stgs->direction</code>.
     * 
     * @param work workspace structure
     * @param i iteration count
     * @return status code; negative status corresponds to error. 
     * 
     * @see direction_type
     */
    scs_int computeDirection(Work *work, scs_int i);


#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

