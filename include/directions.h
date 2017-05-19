#ifndef DIRECTIONS_H
#define DIRECTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"



    /**
     * Adds a new pair (Yk, Sk) into the YS-cache.
     * 
     * @param w
     * @param broydenMemory
     */
    void pushToYSCache(Work *w, YSCache * broydenMemory);

    /**
     * 
     * @param w
     * @param iter
     */
    void computeLBroyden(Work *w, scs_int iter);




#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

