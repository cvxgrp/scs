#include "directions.h"

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