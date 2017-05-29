#include "directions.h"

scs_int resetSUCache(SUCache * cache) {
    cache->mem_cursor = 0; /* set active memory to 0 */
    return SU_CACHE_RESET;
}

scs_int computeLSBroyden(Work *work) {
    /* --- DECLARATIONS --- */
    SUCache * cache; /* the SU cache (pointer) */
    scs_int i; /* index */
    scs_float * s_tilde_current; /* s_tilde (which is updated) */
    scs_float * u_new; /* new value of u */
    scs_float * s_i; /* pointer to the current value of s_i */
    scs_float * u_i; /* pointer to the current value of u_i */
    scs_float ip; /* temporary float to store inner products */
    scs_float gamma; /* scalar gamma as in (6.5e) */
    scs_float theta = 0; /* theta */
    const scs_int l = work->l; /* size of vectors */
    const scs_float theta_bar = work->stgs->thetabar; /* parameter in Powell's trick */

    cache = work->su_cache;

    /* d [work->dir] = -R [work->R] */
    setAsScaledArray(work->dir, work->R, -1.0, l);

    /* s_tilde_current = y [work->Yk] */
    /* use the end of the cache to store s_tilde_current */
    /* later we use the same position of the S-buffer to store the current Sk */
    s_tilde_current = cache->S + (cache->mem_cursor * l);
    memcpy(s_tilde_current, work->Yk, l * sizeof (scs_float));

    /* update s and d */
    for (i = 0; i < cache->mem_cursor; ++i) {
        s_i = cache->S + i * l; /* retrieve s_i from the cache */
        u_i = cache->U + i * l; /* retrieve u_i from the cache */
        ip = innerProd(s_i, s_tilde_current, l);
        addScaledArray(s_tilde_current, u_i, ip, l); /* update s_tilde */
        ip = innerProd(s_i, work->dir, l);
        addScaledArray(work->dir, u_i, ip, l); /* update direction */
    }

    /* compute theta */
    gamma = innerProd(s_tilde_current, work->Sk, l);
    gamma /= calcNormSq(work->Sk, l);

    if (scs_abs(gamma) >= theta_bar) {
        theta = 1;
    } else {
        theta = (1 - scs_sgn(gamma) * theta_bar) / (1 - gamma);
    }


    /* FINALISE */

    /* s_tilde_current = (1-theta)*s + theta*s_tilde_current */
    for (i = 0; i < l; ++i) {
        s_tilde_current[i] = (1 - theta) * work->Sk[i] + theta * s_tilde_current[i];
    }
    /* update u_new (at the end of the buffer) */
    u_new = cache->U + (cache->mem_cursor * l);
    ip = innerProd(work->Sk, s_tilde_current, l);
    for (i = 0; i < l; ++i) {
        u_new[i] = (work->Sk[i] - s_tilde_current[i]) / ip;
    }
    /* update direction */
    ip = innerProd(work->Sk, work->dir, l); /* s'd */
    for (i = 0; i < l; ++i) {
        work->dir[i] += ip * u_new[i];
    }

    /* push s into the buffer */
    memcpy(s_tilde_current, work->Sk, l * sizeof (scs_float));

    cache->mem_cursor++;

    if (cache->mem_cursor >= cache->mem) {
        return resetSUCache(cache); /* returns SU_CACHE_RESET */
    }

    return SU_CACHE_INCREMENT;
}