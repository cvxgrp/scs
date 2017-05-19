#include "directions.h"

void pushToYSCache(Work *w, YSCache * broydenMemory) {
    broydenMemory -> mem_idx++;
    
    if (broydenMemory -> mem_idx >= w->stgs->memory) {
        broydenMemory -> mem_idx = 0;
    }
    broydenMemory->mem_current++; // increase current memory
    if (broydenMemory->mem_current > w->stgs->memory) {
        broydenMemory->mem_current = w->stgs->memory; // buffer is now full
    }

    scs_float *Y = broydenMemory->Y;
    scs_float *yk = w->Yk;
    Y += (broydenMemory -> mem_idx) * w->l;
    memcpy(Y, yk, w->l * sizeof (scs_float));
    
    scs_float *S = broydenMemory->S;
    scs_float *sk = w->Sk;
    S += (broydenMemory -> mem_idx) * w->l;
    memcpy(S, sk, w->l * sizeof (scs_float));
    
    scs_float *YS = broydenMemory->YS;
    YS[broydenMemory -> mem_idx] = innerProd(w->Yk, w->Sk, w->l);
    
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