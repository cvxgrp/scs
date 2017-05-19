#include "directions.h"

/* Broyden for theta=1 */
void computeBroyden(Work *w, scs_int iter) {
    DEBUG_FUNC
    scs_float Ysk;
    scs_int l = w->m + w->n + 1;
    if (iter == 0) {      
        
    }

    if ((w->stgs->sc_init && iter == 1) || w->stgs->tRule == 1
            || w->stgs->tRule == 2) {
        Ysk = innerProd(w->Yk, w->Sk, l);
    }
    if (w->stgs->sc_init && iter == 1 ) {
        
    }
}

void computeLBroyden(Work *w, scs_int iter) {
    DEBUG_FUNC
    scs_float Ysk, qf, theta;
    scs_int l = w->m + w->n + 1;
    scs_int skip = 0;

    Ysk = innerProd(w->Yk, w->Sk, l);

    switch (w->stgs->tRule) {
        case 1:
        case 2:
        case 3:
            qf = -w->u[l - 1] * innerProd(w->Sk, w->sc_R_prev, l);
            if (Ysk < w->stgs->delta * ABS(qf)) {
                theta = (1.0 - SGN(qf) * w->stgs->delta) * qf / (qf - Ysk);
                scaleArray(w->Yk, theta, l);
                addScaledArray(w->Yk, w->sc_R_prev, l, -w->u[l - 1]*(1.0 * theta));
            }
            break;
        case 4:
            if (w->nrmR_con < 1.0)
                w->stgs->alphaC = 3.0;
            if (Ysk / innerProd(w->Sk, w->Sk, l) <=
                    (1e-6) * POWF(w->nrmR_con, w->stgs->alphaC)) {
                skip = 1;
            }
            break;
    }
    if (skip != 0) {
        if (iter < w->stgs->memory) {

        }
    }
}

/* shift columns of a 2D array to the lest by 1 place */
void shiftColLeft(scs_float *arr, scs_int rows, scs_int cols) {
    scs_float temp = arr[0];
    memmove(&arr[0], &arr[1], (rows * cols - 1) * sizeof *arr);
    arr[rows * cols - 1] = temp;
}


// for theta = 1
void computeLMBroyden(Work *w, scs_int iter) {
    
    
}