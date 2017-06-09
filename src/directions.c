#include "directions.h"

static scs_float * HY; /* Vector H*Y of length l */

scs_int resetSUCache(SUCache * cache) {
    cache->mem_cursor = 0; /* set active memory to 0 */
    RETURN SU_CACHE_RESET;
}

scs_int computeLSBroyden(Work *work) {
    /* --- DECLARATIONS --- */
    SUCache * cache; /* the SU cache (pointer) */
    scs_int i; /* index */
    scs_float * s_tilde_current; /* s_tilde (which is updated) */
    scs_float * u_new; /* new value of u */
    scs_float ip; /* temporary float to store inner products */
    scs_float s_norm_sq; /* scalar gamma as in (6.5e) */
    scs_float theta = 0; /* theta */
    const scs_int l = work->l; /* size of vectors */
    const scs_float theta_bar = work->stgs->thetabar; /* parameter in Powell's trick */

    cache = work->su_cache; /* cache of Sk and Uk */

    /* d [work->dir] = -R [work->R] */
    setAsScaledArray(work->dir, work->R, -1.0, l);

    /* s_tilde_current = y [work->Yk]                                           */
    /* use the end of the cache to store s_tilde_current                        */
    /* later we use the same position of the S-buffer to store the current Sk   */
    s_tilde_current = cache->S + (cache->mem_cursor * l);
    memcpy(s_tilde_current, work->Yk, l * sizeof (scs_float));

    /* update s and d */
    for (i = 0; i < cache->mem_cursor; ++i) {
        scs_float * s_i; /* pointer to the current value of s_i */
        scs_float * u_i; /* pointer to the current value of u_i */
        s_i = cache->S + i * l; /* retrieve s_i from the cache */
        u_i = cache->U + i * l; /* retrieve u_i from the cache */
        ip = innerProd(s_i, s_tilde_current, l);
        addScaledArray(s_tilde_current, u_i, l, ip); /* update s_tilde */
        ip = innerProd(s_i, work->dir, l);
        addScaledArray(work->dir, u_i, l, ip); /* update direction */
    }

    /* compute theta */
    ip = innerProd(s_tilde_current, work->Sk, l);
    s_norm_sq = calcNormSq(work->Sk, l);

    if (ABS(ip) >= theta_bar * s_norm_sq) {
        theta = 1;
    } else {
        theta = s_norm_sq * (1 - SGN(ip) * theta_bar) / (s_norm_sq - ip);
        /* s_tilde_current = (1-theta)*s + theta*s_tilde_current */
        for (i = 0; i < l; ++i) {
            s_tilde_current[i] = (1 - theta) * work->Sk[i] + theta * s_tilde_current[i];
        }
    }

    /* FINALISE */

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

    cache->mem_cursor++; /* move the cursor */

    /* if the cursor has exceeded the last position, reset the cache */
    if (cache->mem_cursor >= cache->mem) {
        RETURN resetSUCache(cache); /* returns SU_CACHE_RESET */
    }

    RETURN SU_CACHE_INCREMENT;
}

scs_int computeFullBroyden(Work *work, scs_int i) {
    scs_float ip = 0;
    scs_float tmp = 0;


    if (i == 0 || HY == SCS_NULL) {
        /* HY is allocated the first time this function is called (that is, for i==0) */
        HY = malloc(work->l * sizeof (*HY));
        if (HY == SCS_NULL){
            scs_printf("ERROR: allocating `HY` in `computeFullBroyden` failure\n");
            RETURN DIRECTION_ERROR;
        }
    }

    if ((work->stgs->broyden_init_scaling && i == 1)
            || (work->stgs->tRule == 1 || work->stgs->tRule == 2)) {
        ip = innerProd(work->Yk, work->Sk, work->l);
    }

    if (work->stgs->broyden_init_scaling && i == 1) {
        scs_int i;
        tmp = ip / calcNorm(work->Yk, work->l);
        for (i = 0; i < work->l; ++i) {
            work->H[i * (work->l + 1)] = tmp;
        }
    }

    return 0;
}

scs_int computeDirection(Work *work, scs_int i) {
    scs_int j;
    scs_int status = DIRECTION_SUCCESS;
    if (work->stgs->direction == fixed_point_residual) {
        for (j = 0; j < work->l; ++j) {
            work->dir[j] = -work->R[j];
        }
        status = DIRECTION_SUCCESS;
    } else if (work->stgs->direction == restarted_broyden) {
        status = computeLSBroyden(work);
    } else if (work->stgs->direction == restarted_broyden_v2) {
        status = DIRECTION_ERROR; /* Not implemented yet */
    } else if (work->stgs->direction == full_broyden) {

    }
    RETURN status;
}

void freeFullBroyden() {
    if (HY != SCS_NULL) {
        scs_free(HY);
    }
}