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
     * iteration (current FPR, values of \f$s_k\f$, \f$y_k\f$ etc). 
     * 
     * Directions are computed according to the following procedure:
     * 
     *  1.  Inputs: \f$y=y_k\f$, \f$r=Rx_k\f$, \f$s=s_k\f$, \f$\bar{\theta}\f$, \f$m\f$ (memory)
     *  2.  Buffer: \f$(\mathbf{s}, \mathbf{u})\f$
     *  3.  \f$d_\star \gets -r\f$
     *  4.  \f$s_\star \gets y\f$
     *  5.  for \f$i=0,\ldots, m-1\f$,
     *        1.   \f$s_\star \gets s_\star + \langle \mathbf{s}_i, s_\star\rangle \mathbf{u}_i\f$
     *        2.   \f$d_\star \gets d_\star + \langle \mathbf{s}_i, d_\star\rangle \mathbf{u}_i\f$
     *  6.  \f$\gamma \gets \langle s, s_\star \rangle/\|s\|^2\f$
     *  7.  \f$\theta \gets \begin{cases}1,&\text{if } |\gamma| > \bar{\theta}\\ \frac{1-\mathrm{sgn}(\gamma)\bar{\theta}}{1-\gamma},&\text{otherwise}\end{cases}\f$
     *  8.  \f$s_\star \gets (1-\theta)s + \theta s_\star\f$
     *  9.  \f$u_\star \gets \frac{s-s_\star}{\langle s, s_\star \rangle}\f$ and push it into the buffer
     *  10. \f$d_\star \gets d_\star + \langle s, d_\star\rangle u_\star\f$
     *  11. Add \f$s\f$ into the buffer or empty/reset it if it is full
     * 
     * @return status code of the method.
     */
    scs_int computeLSBroyden(Work *work);


#ifdef __cplusplus
}
#endif

#endif /* DIRECTIONS_H */

