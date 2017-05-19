#include "directions.h"

scs_int pushToYSCache(Work *w, YSCache * ys_cache) {
    
    scs_int cache_status = YS_CACHE_INCREMENT;
    
    /* increment the memory */
    ys_cache -> mem_idx++;        
    if (ys_cache -> mem_idx >= w->stgs->memory) {
        cache_status = YS_CACHE_RESET;
        ys_cache -> mem_idx = 0;
    }
    
    /* update the current memory */
    ys_cache->mem_current++; // increase current memory
    if (ys_cache->mem_current > w->stgs->memory) {
        cache_status = YS_CACHE_FULL;
        ys_cache->mem_current = w->stgs->memory; /* buffer is now full */
    }

    /* push Yk into the right position in the cache */
    scs_float *Y = ys_cache->Y;
    scs_float *yk = w->Yk;
    Y += (ys_cache -> mem_idx) * w->l;
    memcpy(Y, yk, w->l * sizeof (scs_float));
    
    /* push Sk into the right position in the cache */
    scs_float *S = ys_cache->S;
    scs_float *sk = w->Sk;
    S += (ys_cache -> mem_idx) * w->l;
    memcpy(S, sk, w->l * sizeof (scs_float));
    
    /* Compute the inner product and store in the right position in the cache */
    scs_float *YS = ys_cache->YS;
    YS[ys_cache -> mem_idx] = innerProd(w->Yk, w->Sk, w->l);
    
    return cache_status;
}

void resetYSCache(YSCache * broydenMemory){
    
}

void computeLBroyden(Work *w, scs_int iter) {
    DEBUG_FUNC
    scs_float Ysk, qf, theta;
    scs_int l = w->m + w->n + 1;
    scs_int skip = 0;



    Ysk = innerProd(w->Yk, w->Sk, l);



    //
    //    switch (w->stgs->tRule) {
    //        case 1:
    //        case 2:
    //        case 3:
    //            qf = -w->u[l - 1] * innerProd(w->Sk, w->sc_R_prev, l);
    //            if (Ysk < w->stgs->delta * ABS(qf)) {
    //                theta = (1.0 - SGN(qf) * w->stgs->delta) * qf / (qf - Ysk);
    //                scaleArray(w->Yk, theta, l);
    //                addScaledArray(w->Yk, w->sc_R_prev, l, -w->u[l - 1]*(1.0 * theta));
    //            }
    //            break;
    //        case 4:
    //            if (w->nrmR_con < 1.0)
    //                w->stgs->alphaC = 3.0;
    //            if (Ysk / innerProd(w->Sk, w->Sk, l) <=
    //                    (1e-6) * POWF(w->nrmR_con, w->stgs->alphaC)) {
    //                skip = 1;
    //            }
    //            break;
    //    }
    //    if (skip != 0) {
    //        if (iter < w->stgs->memory) {
    //
    //        }
    //    }
}


void computeLMBroyden(Work *w, scs_int iter) {


}