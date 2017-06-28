#include "scs.h"
#include "normalize.h"
#include "directions.h"

#ifndef EXTRAVERBOSE
/* if verbose print summary output every this num iterations */
#define PRINT_INTERVAL 10
/* check for convergence every this num iterations */
#define CONVERGED_INTERVAL 1
#else
#define PRINT_INTERVAL 1
#define CONVERGED_INTERVAL 1
#endif

/* tolerance at which we declare problem indeterminate */
#define INDETERMINATE_TOL 1e-9

timer globalTimer;

/* printing header */
static const char *HEADER[] = {
    " Iter ", " pri res ", " dua res ", " rel gap ",
    " pri obj ", " dua obj ", " kap/tau ", "   FPR   ", " time (s)",
};
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 9;
static const scs_int LINE_LEN = 85;

static scs_int scs_isnan(scs_float x) {
    DEBUG_FUNC
    RETURN(x == NAN || x != x);
}

static SUCache * initSUCache(scs_int memory, scs_int l) {
    SUCache * cache = scs_calloc(1, sizeof (*cache));
    if (cache == SCS_NULL) {
        /* LCOV_EXCL_START */
        scs_printf("ERROR: allocating YSCache failure\n");
        RETURN SCS_NULL;
        /* LCOV_EXCL_STOP */
    }

    /* we allocate one extra memory position because it's needed */
    cache->S = scs_malloc((1 + memory) * l * sizeof (scs_float)); /* S: l-by-mem */
    cache->U = scs_malloc((1 + memory) * l * sizeof (scs_float)); /* U: l-by-mem */


    /* the cache must know its memory length */
    cache->mem = memory;

    /* initial active memory is 0 */
    resetSUCache(cache);
    RETURN cache;
}

static void freeYSCache(SUCache * cache) {
    if (cache == SCS_NULL) {
        return;
    }
    if (cache->S) {
        scs_free(cache->S);
    }
    if (cache->U) {
        scs_free(cache->U);
    }
    scs_free(cache);
    RETURN;
}

static void freeWork(Work *w) {
    DEBUG_FUNC
    if (w == SCS_NULL)
        RETURN;
    if (w->u != SCS_NULL)
        scs_free(w->u);
    if (w->v != SCS_NULL)
        scs_free(w->v);
    if (w->u_t != SCS_NULL)
        scs_free(w->u_t);
    if (w->u_prev != SCS_NULL)
        scs_free(w->u_prev);
    if (w->h != SCS_NULL)
        scs_free(w->h);
    if (w->g != SCS_NULL)
        scs_free(w->g);
    if (w->b != SCS_NULL)
        scs_free(w->b);
    if (w->c != SCS_NULL)
        scs_free(w->c);
    if (w->pr != SCS_NULL)
        scs_free(w->pr);
    if (w->dr != SCS_NULL)
        scs_free(w->dr);
    if (w->scal != SCS_NULL) {
        if (w->scal->D != SCS_NULL)
            scs_free(w->scal->D);
        if (w->scal->E != SCS_NULL)
            scs_free(w->scal->E);
        scs_free(w->scal);
    }
    if (w->u_b != SCS_NULL)
        scs_free(w->u_b);

    if (w->stgs->do_super_scs == 1) {
        if (w->R != SCS_NULL)
            scs_free(w->R);
        if (w->R_prev != SCS_NULL)
            scs_free(w->R_prev);
        if (w->dir != SCS_NULL)
            scs_free(w->dir);
        if (w->dut != SCS_NULL)
            scs_free(w->dut);
        if (w->Sk != SCS_NULL)
            scs_free(w->Sk);
        if (w->Yk != SCS_NULL)
            scs_free(w->Yk);
        if (w->wu != SCS_NULL)
            scs_free(w->wu);
        if (w->wu_t != SCS_NULL)
            scs_free(w->wu_t);
        if (w->wu_b != SCS_NULL)
            scs_free(w->wu_b);
        if (w->Rwu != SCS_NULL)
            scs_free(w->Rwu);
        if (w->su_cache != SCS_NULL)
            freeYSCache(w->su_cache);
        if (w->s_b != SCS_NULL)
            scs_free(w->s_b);
        if (w->H != SCS_NULL) {
            scs_free(w->H);
        }
    }
    scs_free(w);
    RETURN;
}

/* LCOV_EXCL_START */
static void printInitHeader(const Data *d, const Cone *k) {
    DEBUG_FUNC
    scs_int i;
    Settings *stgs = d->stgs;
    char *coneStr = getConeHeader(k);
    char *linSysMethod = getLinSysMethod(d->A, d->stgs);
    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("-");
    }
    scs_printf("\n\tSCS v%s - Splitting Conic Solver\n\t(c) Brendan "
            "O'Donoghue, Stanford University, 2012-2016\n",
            scs_version());
    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("-");
    }
    scs_printf("\n");
    if (linSysMethod) {
        scs_printf("Lin-sys: %s\n", linSysMethod);
        scs_free(linSysMethod);
    }
    if (stgs->normalize) {
        scs_printf("eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i, "
                "scale = %2.2f\n",
                stgs->eps, stgs->alpha, (int) stgs->max_iters,
                (int) stgs->normalize, stgs->scale);
    } else {
        scs_printf("eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i\n",
                stgs->eps, stgs->alpha, (int) stgs->max_iters,
                (int) stgs->normalize);
    }
    scs_printf("Variables n = %i, constraints m = %i\n", (int) d->n, (int) d->m);
    scs_printf("%s", coneStr);
    scs_free(coneStr);
#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
    RETURN;
}

/* LCOV_EXCL_STOP */

static void populateOnFailure(scs_int m, scs_int n, Sol *sol, Info *info,
        scs_int statusVal, const char *msg) {
    DEBUG_FUNC
    if (info) {
        info->relGap = NAN;
        info->resPri = NAN;
        info->resDual = NAN;
        info->pobj = NAN;
        info->dobj = NAN;
        info->iter = -1;
        info->statusVal = statusVal;
        info->solveTime = NAN;
        strcpy(info->status, msg);
    }
    if (sol) {
        if (n > 0) {
            if (sol->x == SCS_NULL)
                sol->x = scs_malloc(sizeof (scs_float) * n);
            scaleArray(sol->x, NAN, n);
        }
        if (m > 0) {
            if (sol->y == SCS_NULL)
                sol->y = scs_malloc(sizeof (scs_float) * m);
            scaleArray(sol->y, NAN, m);
            if (sol->s == SCS_NULL)
                sol->s = scs_malloc(sizeof (scs_float) * m);
            scaleArray(sol->s, NAN, m);
        }
    }
    RETURN;
}

static scs_int failure(Work *w, scs_int m, scs_int n, Sol *sol, Info *info,
        scs_int stint, const char *msg, const char *ststr) {
    DEBUG_FUNC
    scs_int status = stint;
    populateOnFailure(m, n, sol, info, status, ststr);
    scs_printf("Failure:%s\n", msg);
    endInterruptListener();
    RETURN status;
}

static void warmStartVars(Work *w, const Sol *sol) {
    DEBUG_FUNC
    scs_int i, n = w->n, m = w->m;
    memset(w->v, 0, n * sizeof (scs_float));
    memcpy(w->u, sol->x, n * sizeof (scs_float));
    memcpy(&(w->u[n]), sol->y, m * sizeof (scs_float));
    if (!w->stgs->do_super_scs) {
        memcpy(&(w->v[n]), sol->s, m * sizeof (scs_float));
        w->v[n + m] = 0.0;
    }
    w->u[n + m] = 1.0;
#ifndef NOVALIDATE
    for (i = 0; i < n + m + 1; ++i) {
        if (scs_isnan(w->u[i])) {
            w->u[i] = 0;
        }
        if (scs_isnan(w->v[i])) {
            w->v[i] = 0;
        }
    }
#endif
    if (w->stgs->normalize) {
        normalizeWarmStart(w);
    }
    RETURN;
}

static scs_float calcPrimalResid(Work *w, const scs_float *x,
        const scs_float *s, const scs_float tau,
        scs_float *nmAxs) {
    DEBUG_FUNC
    scs_int i;
    scs_float pres = 0, scale, *pr = w->pr;
    *nmAxs = 0;
    memset(pr, 0, w->m * sizeof (scs_float));
    accumByA(w->A, w->p, x, pr);
    addScaledArray(pr, s, w->m, 1.0); /* pr = Ax + s */
    for (i = 0; i < w->m; ++i) {
        scale =
                w->stgs->normalize ? w->scal->D[i] / (w->sc_b * w->stgs->scale) : 1;
        scale = scale * scale;
        *nmAxs += (pr[i] * pr[i]) * scale;
        pres += (pr[i] - w->b[i] * tau) * (pr[i] - w->b[i] * tau) * scale;
    }
    *nmAxs = SQRTF(*nmAxs);
    RETURN SQRTF(pres); /* norm(Ax + s - b * tau), for superSCS: norm(Ax_b + s_b - b * tau_b) */
}

static scs_float calcDualResid(Work *w, const scs_float *y, const scs_float tau,
        scs_float *nmATy) {
    DEBUG_FUNC
    scs_int i;
    scs_float dres = 0, scale, *dr = w->dr;
    *nmATy = 0;
    memset(dr, 0, w->n * sizeof (scs_float));
    accumByAtrans(w->A, w->p, y, dr); /* dr = A'y */
    for (i = 0; i < w->n; ++i) {
        scale =
                w->stgs->normalize ? w->scal->E[i] / (w->sc_c * w->stgs->scale) : 1;
        scale = scale * scale;
        *nmATy += (dr[i] * dr[i]) * scale;
        dres += (dr[i] + w->c[i] * tau) * (dr[i] + w->c[i] * tau) * scale;
    }
    *nmATy = SQRTF(*nmATy);
    RETURN SQRTF(dres); /* norm(A'y + c * tau), for superSCS: norm(A'y_b + c * tau_b)*/
}

/* calculates un-normalized quantities */
static void calcResiduals(Work *w, struct residuals *r, scs_int iter) {
    DEBUG_FUNC
    scs_float *x;
    scs_float *y;
    scs_float *s;
    scs_float nmpr_tau;
    scs_float nmdr_tau;
    scs_float nmAxs_tau;
    scs_float nmATy_tau;
    scs_float cTx, bTy;
    scs_int n = w->n, m = w->m;

    /* checks if the residuals are unchanged by checking iteration */
    if (r->lastIter == iter) {
        RETURN;
    }
    r->lastIter = iter;


    if (!w->stgs->do_super_scs) {
        s = &(w->v[w->n]);
        x = w->u;
        y = &(w->u[w->n]);

        r->tau = ABS(w->u[n + m]);
        r->kap = ABS(w->v[n + m]) /
                (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);
    } else {
        s = w->s_b;
        x = w->u_b;
        y = &(w->u_b[w->n]);

        r->kap = w->kap_b;
        r->tau = w->u_b[n + m]; /* it's actually tau_b */
    }

    nmpr_tau = calcPrimalResid(w, x, s, r->tau, &nmAxs_tau);
    nmdr_tau = calcDualResid(w, y, r->tau, &nmATy_tau);

    r->bTy_by_tau =
            innerProd(y, w->b, m) /
            (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);
    r->cTx_by_tau =
            innerProd(x, w->c, n) /
            (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);

    r->resInfeas =
            r->bTy_by_tau < 0 ? w->nm_b * nmATy_tau / -r->bTy_by_tau : NAN;
    r->resUnbdd =
            r->cTx_by_tau < 0 ? w->nm_c * nmAxs_tau / -r->cTx_by_tau : NAN;

    bTy = r->bTy_by_tau / r->tau;
    cTx = r->cTx_by_tau / r->tau;

    r->resPri = nmpr_tau / (1 + w->nm_b) / r->tau;
    r->resDual = nmdr_tau / (1 + w->nm_c) / r->tau;
    r->relGap = ABS(cTx + bTy) / (1 + ABS(cTx) + ABS(bTy));
    RETURN;
}

static void calcResidualsSuperscs(
        Work *w,
        struct residuals *r,
        scs_int iter) {
    DEBUG_FUNC
    scs_float *xb;
    scs_float *yb;
    scs_float *sb;
    scs_float cTx;
    scs_float bTy;
    scs_float *pr = w->pr;
    scs_float *dr = w->dr;
    scs_int n = w->n;
    scs_int m = w->m;
    scs_int i;
    scs_float norm_D_Axs; /* norm of D*(Ax+s), intermediate variable */
    scs_float norm_E_ATy; /* norm of E*A'*y,   intermediate variable */
    scs_float tmp_cTx; /* c'x */
    scs_float tmp_bTy; /* b'y */
    const scs_float temp1 = w->sc_b * w->stgs->scale; /* auxiliary variable #1 */
    const scs_float temp2 = w->sc_c * temp1; /* auxiliary variable #2 */
    const scs_float temp3 = w->sc_c * w->stgs->scale; /* auxiliary variable #3 */

    /* checks if the residuals are unchanged by checking iteration */
    if (r->lastIter == iter) {
        RETURN;
    }
    r->lastIter = iter;

    sb = w->s_b;
    xb = w->u_b;
    yb = &(w->u_b[n]);

    r->kap = w->kap_b;
    r->tau = w->u_b[n + m]; /* it's actually tau_b */

    memset(pr, 0, w->m * sizeof (scs_float)); /* pr = 0 */
    memset(dr, 0, w->n * sizeof (scs_float)); /* dr = 0 */

    accumByA(w->A, w->p, xb, pr); /* pr = A xb */
    addScaledArray(pr, sb, w->m, 1.0); /* pr = A xb + sb */
    /* --- compute ||D(Ax + s)|| --- */
    norm_D_Axs = 0;
    for (i = 0; i < m; ++i) {
        scs_float tmp = w->scal->D[i] * pr[i];
        norm_D_Axs += tmp;
    }
    norm_D_Axs = SQRTF(norm_D_Axs);
    addScaledArray(pr, w->b, m, -r->tau); /* pr = A xb + sb - b taub */

    accumByAtrans(w->A, w->p, yb, dr); /* dr = A' yb */
    /* --- compute ||E A' yb|| --- */
    norm_E_ATy = 0;
    for (i = 0; i < n; ++i) {
        scs_float tmp = w->scal->E[i] * dr[i];
        norm_E_ATy += tmp;
    }
    norm_E_ATy = SQRTF(norm_E_ATy);
    addScaledArray(dr, w->c, w->n, r->tau); /* dr = A' yb + c taub */

    /*
     * bTy_by_tau = b'yb / (scale*sc_c*sc_b)
     * cTx_by_tau = c'xb / (scale*sc_c*sc_b)
     */
    tmp_bTy = innerProd(yb, w->b, m);
    r->bTy_by_tau = tmp_bTy / (w->stgs->normalize ? (temp2) : 1);
    tmp_cTx = innerProd(xb, w->c, n);
    r->cTx_by_tau = tmp_cTx / (w->stgs->normalize ? (temp2) : 1);

    /*
     * bTy = b'yb / (scale*sc_c*sc_b) / taub
     * cTx = c'xb / (scale*sc_c*sc_b) / taub
     */
    bTy = r->bTy_by_tau / r->tau;
    cTx = r->cTx_by_tau / r->tau;

    /* PRIMAL RESIDUAL */
    if (w->stgs->normalize) {
        r->resPri = 0;
        for (i = 0; i < m; ++i) {
            scs_float tmp = w->scal->D[i] * pr[i];
            r->resPri += tmp * tmp;
        }
        r->resPri = SQRTF(r->resPri) / r->tau;
        r->resPri /= ((1 + w->nm_b) * temp1);
    } else {
        r->resPri = calcNorm(pr, m) / r->tau;
        r->resPri /= (1 + w->nm_b);
    }

    /* DUAL RESIDUAL */
    if (w->stgs->normalize) {
        r->resDual = 0;
        for (i = 0; i < n; ++i) {
            scs_float tmp = w->scal->E[i] * dr[i];
            r->resDual += tmp * tmp;
        }
        r->resDual = SQRTF(r->resDual) / r->tau;
        r->resDual /= ((1 + w->nm_c) * temp3);
    } else {
        r->resDual = calcNorm(dr, n) / r->tau;
        r->resDual /= (1 + w->nm_c);
    }

    /* UNBOUNDEDNESS */
    if (tmp_cTx < 0) {
        scs_float norm_Ec = 0;
        for (i = 0; i < n; ++i) {
            scs_float tmp = w->scal->E[i] * w->c[i];
            norm_Ec += tmp * tmp;
        }
        r->resUnbdd = -SQRTF(norm_Ec) * norm_D_Axs / tmp_cTx;
        r->resUnbdd /= w->stgs->normalize ? w->stgs->scale : 1;
    } else {
        r->resUnbdd = NAN;
    }

    /* INFEASIBILITY */
    if (tmp_bTy < 0) {
        scs_float norm_Db = 0;
        for (i = 0; i < m; ++i) {
            scs_float tmp = w->scal->D[i] * w->b[i];
            norm_Db += tmp * tmp;
        }
        r->resInfeas = -SQRTF(norm_Db) * norm_E_ATy / tmp_bTy;
        r->resInfeas /= w->stgs->normalize ? w->stgs->scale : 1;
    } else {
        r->resInfeas = NAN;
    }

    r->relGap = ABS(cTx + bTy) / (1 + ABS(cTx) + ABS(bTy));
    RETURN;
}

static void coldStartVars(Work *w) {
    DEBUG_FUNC
    scs_int l = w->l;
    memset(w->u, 0, l * sizeof (scs_float));
    w->u[l - 1] = SQRTF((scs_float) l);
    if (!w->stgs->do_super_scs) {
        memset(w->v, 0, l * sizeof (scs_float));
        w->v[l - 1] = SQRTF((scs_float) l);
    }
    RETURN;
}

/* status < 0 indicates failure */
static scs_int projectLinSys(Work *w, scs_int iter) {
    /* ut = u + v */
    DEBUG_FUNC
    scs_int n = w->n, m = w->m, l = n + m + 1, status;
    memcpy(w->u_t, w->u, l * sizeof (scs_float));
    addScaledArray(w->u_t, w->v, l, 1.0);

    scaleArray(w->u_t, w->stgs->rho_x, n);

    addScaledArray(w->u_t, w->h, l - 1, -w->u_t[l - 1]);
    addScaledArray(w->u_t, w->h, l - 1,
            -innerProd(w->u_t, w->g, l - 1) / (w->gTh + 1));
    scaleArray(&(w->u_t[n]), -1, m);

    status = solveLinSys(w->A, w->stgs, w->p, w->u_t, w->u, iter);

    w->u_t[l - 1] += innerProd(w->u_t, w->h, l - 1);

    RETURN status;
}

/* status < 0 indicates failure */
static scs_int projectLinSysv2(scs_float * u_t, scs_float * u, Work * w, scs_int iter) {
    DEBUG_FUNC;
    const scs_int l = w->l;
    scs_int status;

    /*TODO make more efficient; eliminate memcpy? use setAsScaledArray? */
    
    /* ut(1:n) = rho_x*ut(1:n); */
    memcpy(u_t + w->n, u + w->n, (w->m + 1) * sizeof (scs_float));
    setAsScaledArray(u_t, u, w->stgs->rho_x, w->n);

    /* ut(1:l-1) = ut(1:l-1) - ut(end) * h      */
    addScaledArray(u_t, w->h, l - 1, -u_t[l - 1]);

    /* ut -= scalar * h                         */
    addScaledArray(u_t, w->h, l - 1,
            -innerProd(u_t, w->g, l - 1) / (w->gTh + 1));

    /* ut(n+1:end-1) = -ut(n+1:end-1);          */
    scaleArray(u_t + w->n, -1, w->m);

    /* call `solveLinSys` to update ut(1:n+m)   */
    status = solveLinSys(w->A, w->stgs, w->p, u_t, u, iter);

    /* ut(end) = (ut(end) + h'*ut(1:l-1))       */
    u_t[l - 1] += innerProd(u_t, w->h, l - 1);

    RETURN status;
}

/* LCOV_EXCL_START */
void printSol(Work *w, Sol *sol, Info *info) {
    DEBUG_FUNC
    scs_int i;
    scs_printf("%s\n", info->status);
    if (sol->x != SCS_NULL) {
        for (i = 0; i < w->n; ++i) {
            scs_printf("x[%i] = %4f\n", (int) i, sol->x[i]);
        }
    }
    if (sol->y != SCS_NULL) {
        for (i = 0; i < w->m; ++i) {
            scs_printf("y[%i] = %4f\n", (int) i, sol->y[i]);
        }
    }
    RETURN;
}

/* LCOV_EXCL_STOP */

static void updateDualVars(Work *w) {
    DEBUG_FUNC
    scs_int i, n = w->n, l = n + w->m + 1;
    /* this does not relax 'x' variable */
    for (i = n; i < l; ++i) {
        w->v[i] += (w->u[i] - w->stgs->alpha * w->u_t[i] -
                (1.0 - w->stgs->alpha) * w->u_prev[i]);
    }
    RETURN;
}

/* Calculates the fixed point residual R */
static void calcFPRes(scs_float *R, scs_float *u_t, scs_float *u_b, scs_int l) {
    DEBUG_FUNC
    scs_int i;
    for (i = 0; i < l; ++i) {
        R[i] = u_t[i] - u_b[i];
    }
    RETURN;
}

/* status < 0 indicates failure */
static scs_int projectCones(Work *w, const Cone *k, scs_int iter) {
    DEBUG_FUNC
    scs_int i, n = w->n, l = n + w->m + 1, status;
    /* this does not relax 'x' variable */
    for (i = 0; i < n; ++i) {
        w->u[i] = w->u_t[i] - w->v[i];
    }
    for (i = n; i < l; ++i) {
        w->u[i] = w->stgs->alpha * w->u_t[i] +
                (1 - w->stgs->alpha) * w->u_prev[i] - w->v[i];
    }
    /* u = [x;y;tau] */
    status = projDualCone(&(w->u[n]), k, w->coneWork, &(w->u_prev[n]), iter);
    if (w->u[l - 1] < 0.0)
        w->u[l - 1] = 0.0;

    RETURN status;
}

/* status < 0 indicates failure */
static scs_int projectConesv2(scs_float *u_b, scs_float *u_t, scs_float *u, Work *w, const Cone *k, scs_int iter) {
    DEBUG_FUNC
    scs_int i, n = w->n, l = n + w->m + 1, status;
    /* this does not relax 'x' variable */
    for (i = 0; i < l; ++i) {
        u_b[i] = 2 * u_t[i] - u[i];
    }

    /* u = [x;y;tau] */
    status = projDualCone(&(u_b[n]), k, w->coneWork, &(w->u_prev[n]), iter);
    if (u_b[l - 1] < 0.0) {
        u_b[l - 1] = 0.0;
    }
    RETURN status;
}

static scs_int indeterminate(Work *w, Sol *sol, Info *info) {
    DEBUG_FUNC
    strcpy(info->status, "Indeterminate");
    scaleArray(sol->x, NAN, w->n);
    scaleArray(sol->y, NAN, w->m);
    scaleArray(sol->s, NAN, w->m);
    RETURN SCS_INDETERMINATE;
}

static scs_int solved(Work *w, Sol *sol, Info *info, scs_float tau) {
    DEBUG_FUNC
    scaleArray(sol->x, 1.0 / tau, w->n);
    scaleArray(sol->y, 1.0 / tau, w->m);
    scaleArray(sol->s, 1.0 / tau, w->m);
    if (info->statusVal == 0) {
        strcpy(info->status, "Solved/Inaccurate");
        RETURN SCS_SOLVED_INACCURATE;
    }
    strcpy(info->status, "Solved");
    RETURN SCS_SOLVED;
}

static scs_int infeasible(Work *w, Sol *sol, Info *info, scs_float bTy) {
    DEBUG_FUNC
    scaleArray(sol->y, -1 / bTy, w->m);
    scaleArray(sol->x, NAN, w->n);
    scaleArray(sol->s, NAN, w->m);
    if (info->statusVal == 0) {
        strcpy(info->status, "Infeasible/Inaccurate");
        RETURN SCS_INFEASIBLE_INACCURATE;
    }
    strcpy(info->status, "Infeasible");
    RETURN SCS_INFEASIBLE;
}

static scs_int unbounded(Work *w, Sol *sol, Info *info, scs_float cTx) {
    DEBUG_FUNC
    scaleArray(sol->x, -1 / cTx, w->n);
    scaleArray(sol->s, -1 / cTx, w->m);
    scaleArray(sol->y, NAN, w->m);
    if (info->statusVal == 0) {
        strcpy(info->status, "Unbounded/Inaccurate");
        RETURN SCS_UNBOUNDED_INACCURATE;
    }
    strcpy(info->status, "Unbounded");
    RETURN SCS_UNBOUNDED;
}

static void sety(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->y) {
        sol->y = scs_malloc(sizeof (scs_float) * w->m);
    }
    if (!w->stgs->do_super_scs) {
        memcpy(sol->y, &(w->u[w->n]), w->m * sizeof (scs_float));
    } else {
        memcpy(sol->y, &(w->u_b[w->n]), w->m * sizeof (scs_float));
    }
    RETURN;
}

static void sets(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->s) {
        sol->s = scs_malloc(sizeof (scs_float) * w->m);
    }
    if (!w->stgs->do_super_scs) {
        memcpy(sol->s, &(w->v[w->n]), w->m * sizeof (scs_float));
    } else {
        memcpy(sol->s, w->s_b, w->m * sizeof (scs_float));
    }
    RETURN;
}

static void setx(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->x)
        sol->x = scs_malloc(sizeof (scs_float) * w->n);
    if (!w->stgs->do_super_scs) {
        memcpy(sol->x, w->u, w->n * sizeof (scs_float));
    } else {
        memcpy(sol->x, w->u_b, w->n * sizeof (scs_float));
    }
    RETURN;
}

scs_int isSolvedStatus(scs_int status) {
    RETURN status == SCS_SOLVED || status == SCS_SOLVED_INACCURATE;
}

scs_int isInfeasibleStatus(scs_int status) {
    RETURN status == SCS_INFEASIBLE || status == SCS_INFEASIBLE_INACCURATE;
}

scs_int isUnboundedStatus(scs_int status) {
    RETURN status == SCS_UNBOUNDED || status == SCS_UNBOUNDED_INACCURATE;
}

static void getInfo(Work *w, Sol *sol, Info *info, struct residuals *r,
        scs_int iter) {
    DEBUG_FUNC
    info->iter = iter;
    info->resInfeas = r->resInfeas;
    info->resUnbdd = r->resUnbdd;
    if (isSolvedStatus(info->statusVal)) {
        info->relGap = r->relGap;
        info->resPri = r->resPri;
        info->resDual = r->resDual;
        info->pobj = r->cTx_by_tau / r->tau;
        info->dobj = -r->bTy_by_tau / r->tau;
    } else if (isUnboundedStatus(info->statusVal)) {
        info->relGap = NAN;
        info->resPri = NAN;
        info->resDual = NAN;
        info->pobj = -INFINITY;
        info->dobj = -INFINITY;
    } else if (isInfeasibleStatus(info->statusVal)) {
        info->relGap = NAN;
        info->resPri = NAN;
        info->resDual = NAN;
        info->pobj = INFINITY;
        info->dobj = INFINITY;
    }
    RETURN;
}

/* sets solutions, re-scales by inner prods if infeasible or unbounded */
static void getSolution(Work *work, Sol *sol, Info *info, struct residuals *r,
        scs_int iter) {
    DEBUG_FUNC
    scs_int l = work->l;
    if (!work->stgs->do_super_scs) {
        calcResiduals(work, r, iter);
    } else {
        calcResidualsSuperscs(work, r, iter);
        r->kap = ABS(work->kap_b) /
                (work->stgs->normalize ? (work->stgs->scale * work->sc_c * work->sc_b) : 1.0);
    }
    setx(work, sol);
    sety(work, sol);
    sets(work, sol);
    if (info->statusVal == SCS_UNFINISHED) {
        /* not yet converged, take best guess */
        if (r->tau > INDETERMINATE_TOL && r->tau > r->kap) {
            info->statusVal = solved(work, sol, info, r->tau);
        } else if (calcNorm(work->u, l) <
                INDETERMINATE_TOL * SQRTF((scs_float) l)) {
            info->statusVal = indeterminate(work, sol, info);
        } else if (r->bTy_by_tau < r->cTx_by_tau) {
            info->statusVal = infeasible(work, sol, info, r->bTy_by_tau);
        } else {
            info->statusVal = unbounded(work, sol, info, r->cTx_by_tau);
        }
    } else if (isSolvedStatus(info->statusVal)) {
        info->statusVal = solved(work, sol, info, r->tau);
    } else if (isInfeasibleStatus(info->statusVal)) {
        info->statusVal = infeasible(work, sol, info, r->bTy_by_tau);
    } else {
        info->statusVal = unbounded(work, sol, info, r->cTx_by_tau);
    }
    if (work->stgs->normalize) {
        unNormalizeSol(work, sol);
    }
    getInfo(work, sol, info, r, iter);
    RETURN;
}

/* LCOV_EXCL_START */
static void printSummary(Work *w, scs_int i, struct residuals *r,
        timer *solveTimer) {
    DEBUG_FUNC
    scs_printf("%*i|", (int) strlen(HEADER[0]), (int) i);
    scs_printf("%*.2e ", (int) HSPACE, r->resPri);
    scs_printf("%*.2e ", (int) HSPACE, r->resDual);
    scs_printf("%*.2e ", (int) HSPACE, r->relGap);
    scs_printf("%*.2e ", (int) HSPACE, r->cTx_by_tau / r->tau);
    scs_printf("%*.2e ", (int) HSPACE, -r->bTy_by_tau / r->tau);
    scs_printf("%*.2e ", (int) HSPACE, r->kap / r->tau);
    if (w->stgs->do_super_scs) {
        scs_printf("%*.2e ", (int) HSPACE, w->nrmR_con);
    }
    scs_printf("%*.2e ", (int) HSPACE, tocq(solveTimer) / 1e3);
    scs_printf("\n");

#if EXTRAVERBOSE > 0
    scs_printf("Norm u = %4f, ", calcNorm(w->u, w->n + w->m + 1));
    scs_printf("Norm u_t = %4f, ", calcNorm(w->u_t, w->n + w->m + 1));
    if (!w->stgs->do_super_scs) {
        scs_printf("Norm v = %4f, ", calcNorm(w->v, w->n + w->m + 1));
    }
    scs_printf("tau = %4f, ", w->u[w->n + w->m]);
    if (!w->stgs->do_super_scs) {
        scs_printf("kappa = %4f, ", w->v[w->n + w->m]);
    }
    scs_printf("|u - u_prev| = %1.2e, ",
            calcNormDiff(w->u, w->u_prev, w->n + w->m + 1));
    scs_printf("|u - u_t| = %1.2e, ",
            calcNormDiff(w->u, w->u_t, w->n + w->m + 1));
    scs_printf("resInfeas = %1.2e, ", r->resInfeas);
    scs_printf("resUnbdd = %1.2e\n", r->resUnbdd);
#endif

#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
    RETURN;
}

static void printHeader(Work *w, const Cone *k) {
    DEBUG_FUNC
    scs_int i;
    if (w->stgs->warm_start)
        scs_printf("SCS using variable warm-starting\n");
    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("-");
    }
    scs_printf("\n");
    for (i = 0; i < HEADER_LEN - 2; ++i) {
        scs_printf("%s|", HEADER[i]);
    }
    if (w->stgs->do_super_scs) {
        scs_printf("%s|", HEADER[HEADER_LEN - 2]);
    }
    scs_printf("%s\n", HEADER[HEADER_LEN - 1]);
    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("-");
    }
    scs_printf("\n");
#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
    RETURN;
}

/* LCOV_EXCL_STOP */

scs_float getDualConeDist(const scs_float *y, const Cone *k, ConeWork *c,
        scs_int m) {
    DEBUG_FUNC
    scs_float dist;
    scs_float *t = scs_malloc(sizeof (scs_float) * m);
    memcpy(t, y, m * sizeof (scs_float));
    projDualCone(t, k, c, SCS_NULL, -1);
    dist = calcNormInfDiff(t, y, m);
#if EXTRAVERBOSE > 0
    printArray(y, m, "y");
    printArray(t, m, "projY");
    scs_printf("dist = %4f\n", dist);
#endif
    scs_free(t);
    RETURN dist;
}

/* via moreau */
scs_float getPriConeDist(const scs_float *s, const Cone *k, ConeWork *c,
        scs_int m) {
    DEBUG_FUNC
    scs_float dist;
    scs_float *t = scs_malloc(sizeof (scs_float) * m);
    memcpy(t, s, m * sizeof (scs_float));
    scaleArray(t, -1.0, m);
    projDualCone(t, k, c, SCS_NULL, -1);
    dist = calcNormInf(t, m); /* ||s - Pi_c(s)|| = ||Pi_c*(-s)|| */
#if EXTRAVERBOSE > 0
    printArray(s, m, "s");
    printArray(t, m, "(s - ProjS)");
    scs_printf("dist = %4f\n", dist);
#endif
    scs_free(t);
    RETURN dist;
}

/* LCOV_EXCL_START */
static void printFooter(const Data *d, const Cone *k, Sol *sol, Work *w,
        Info *info) {
    DEBUG_FUNC
    scs_int i;
    char *linSysStr = getLinSysSummary(w->p, info);
    char *coneStr = getConeSummary(info, w->coneWork);
    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("-");
    }
    scs_printf("\nStatus: %s\n", info->status);
    if (info->iter == w->stgs->max_iters) {
        scs_printf("Hit max_iters, solution may be inaccurate\n");
    }
    scs_printf("Timing: Solve time: %1.2es\n", info->solveTime / 1e3);

    if (linSysStr) {
        scs_printf("%s", linSysStr);
        scs_free(linSysStr);
    }

    if (coneStr) {
        scs_printf("%s", coneStr);
        scs_free(coneStr);
    }

    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("-");
    }
    scs_printf("\n");

    if (isInfeasibleStatus(info->statusVal)) {
        scs_printf("Certificate of primal infeasibility:\n");
        scs_printf("dist(y, K*) = %.4e\n",
                getDualConeDist(sol->y, k, w->coneWork, d->m));
        scs_printf("|A'y|_2 * |b|_2 = %.4e\n", info->resInfeas);
        scs_printf("b'y = %.4f\n", innerProd(d->b, sol->y, d->m));
    } else if (isUnboundedStatus(info->statusVal)) {
        scs_printf("Certificate of dual infeasibility:\n");
        scs_printf("dist(s, K) = %.4e\n",
                getPriConeDist(sol->s, k, w->coneWork, d->m));
        scs_printf("|Ax + s|_2 * |c|_2 = %.4e\n", info->resUnbdd);
        scs_printf("c'x = %.4f\n", innerProd(d->c, sol->x, d->n));
    } else {
        scs_printf("Error metrics:\n");
        scs_printf("dist(s, K) = %.4e, dist(y, K*) = %.4e, s'y/|s||y| = %.4e\n",
                getPriConeDist(sol->s, k, w->coneWork, d->m),
                getDualConeDist(sol->y, k, w->coneWork, d->m),
                innerProd(sol->s, sol->y, d->m) / calcNorm(sol->s, d->m) /
                calcNorm(sol->y, d->m));
        scs_printf("|Ax + s - b|_2 / (1 + |b|_2) = %.4e\n", info->resPri);
        scs_printf("|A'y + c|_2 / (1 + |c|_2) = %.4e\n", info->resDual);
        scs_printf("|c'x + b'y| / (1 + |c'x| + |b'y|) = %.4e\n", info->relGap);
        for (i = 0; i < LINE_LEN; ++i) {
            scs_printf("-");
        }
        scs_printf("\n");
        scs_printf("c'x = %.4f, -b'y = %.4f\n", info->pobj, info->dobj);
    }
    for (i = 0; i < LINE_LEN; ++i) {
        scs_printf("=");
    }
    scs_printf("\n");
#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
    RETURN;
}

/* LCOV_EXCL_STOP */

static scs_int hasConverged(Work *w, struct residuals *r, scs_int iter) {
    DEBUG_FUNC
    scs_float eps = w->stgs->eps;
    if (r->resPri < eps && r->resDual < eps && r->relGap < eps) {
        RETURN SCS_SOLVED;
    }
    if (r->resUnbdd < eps) {
        RETURN SCS_UNBOUNDED;
    }
    if (r->resInfeas < eps) {
        RETURN SCS_INFEASIBLE;
    }
    RETURN 0;
}

static scs_int validate(const Data *d, const Cone *k) {
    DEBUG_FUNC
    Settings *stgs = d->stgs;
    if (d->m <= 0 || d->n <= 0) {
        /* LCOV_EXCL_START */
        scs_printf("m and n must both be greater than 0; m = %li, n = %li\n",
                (long) d->m, (long) d->n);
        RETURN - 1;
        /* LCOV_EXCL_STOP */
    }
    if (d->m < d->n) {
        /* LCOV_EXCL_START */
        scs_printf("WARN: m less than n, problem likely degenerate\n");
        /* RETURN -1; */
        /* LCOV_EXCL_STOP */
    }
    if (validateLinSys(d->A) < 0) {
        /* LCOV_EXCL_START */
        scs_printf("invalid linear system input data\n");
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (validateCones(d, k) < 0) {
        /* LCOV_EXCL_START */
        scs_printf("cone validation error\n");
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (stgs->max_iters <= 0) {
        /* LCOV_EXCL_START */
        scs_printf("max_iters must be positive (max_iters=%d)\n", stgs->max_iters);
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (stgs->eps <= 0) {
        /* LCOV_EXCL_START */
        scs_printf("eps tolerance must be positive (eps=%g)\n", stgs->eps);
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (stgs->alpha <= 0 || stgs->alpha >= 2) {
        /* LCOV_EXCL_START */
        scs_printf("alpha must be in (0,2) (alpha=%g)\n", stgs->alpha);
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (stgs->rho_x <= 0) {
        /* LCOV_EXCL_START */
        scs_printf("rho_x must be positive (1e-3 works well) (rho_x=%g).\n", stgs->rho_x);
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (stgs->scale <= 0) {
        /* LCOV_EXCL_START */
        scs_printf("Parameter `scale` must be positive (1 works well).\n");
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    if (stgs->do_super_scs != 0 && stgs->do_super_scs != 1) {
        /* LCOV_EXCL_START */
        scs_printf("do_super_scs (=%d) can be either 0 or 1.\n", stgs->do_super_scs);
        RETURN SCS_FAILED;
        /* LCOV_EXCL_STOP */
    }
    /* validate settings related to SuperSCS */
    if (stgs->do_super_scs == 1) {
        if (stgs->thetabar < 0 || stgs->thetabar > 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameters `thetabar` must be a scalar between 0 and 1 (thetabar=%g)\n", stgs->thetabar);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if ((stgs->direction == restarted_broyden || stgs->direction == restarted_broyden_v2)
                && stgs->memory <= 1) {
            /* LCOV_EXCL_START */
            scs_printf("Quasi-Newton memory length (mem=%d) is too low; choose an integer at least equal to 2.\n", stgs->memory);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->beta >= 1 || stgs->beta <= 0) {
            /* LCOV_EXCL_START */
            scs_printf("Stepsize reduction factor (beta=%g) out of bounds.\n", stgs->beta);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->ls < 0) {
            /* LCOV_EXCL_START */
            scs_printf("Illegal maximum number of line search iterations (ls=%d).\n", stgs->ls);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->sigma < 0) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter sigma of the line search (sigma=%g) cannot be negative.\n", stgs->sigma);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->c_bl < 0 || stgs->c_bl >= 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter (c_0=%g) for blind updates out of bounds.\n", stgs->c_bl);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->c1 < 0 || stgs->c1 >= 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter (c1=%g) for step K1 out of bounds.\n", stgs->c1);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->sse < 0 || stgs->sse >= 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter (sse=%g) for step K1 out of bounds.\n", stgs->sse);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->k0 != 0 && stgs->k0 != 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter (k0=%d) can be eiter 0 (k0: off) or 1 (k0: on).\n", stgs->k0);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->k1 != 0 && stgs->k1 != 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter (k1=%d) can be eiter 0 (k1: off) or 1 (k1: on).\n", stgs->k0);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
        if (stgs->k2 != 0 && stgs->k2 != 1) {
            /* LCOV_EXCL_START */
            scs_printf("Parameter (k2=%d) can be eiter 0 (k2: off) or 1 (k2: on).\n", stgs->k0);
            RETURN SCS_FAILED;
            /* LCOV_EXCL_STOP */
        }
    }
    RETURN 0;
}

static Work *initWork(const Data *d, const Cone *k) {
    DEBUG_FUNC
    Work *w = scs_calloc(1, sizeof (*w));
    const scs_int l = d->n + d->m + 1;
    if (d->stgs->verbose) {
        printInitHeader(d, k);
    }
    if (!w) {
        scs_printf("ERROR: allocating work failure\n");
        RETURN SCS_NULL;
    }

    /* get settings and dims from data struct */
    w->stgs = d->stgs;
    w->m = d->m;
    w->n = d->n;
    w->l = l; /* total dimension */

    /* -------------------------------------
     * Workspace allocation:
     * 
     * After every call to scs_malloc or scs_calloc
     * we check whether the allocation has been
     * successful.
     * ------------------------------------- */
    w->u = scs_calloc(l, sizeof (scs_float));
    if (w->u == SCS_NULL) {
        scs_printf("ERROR: `u` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->u_b = scs_calloc(l, sizeof (scs_float));
    if (w->u_b == SCS_NULL) {
        scs_printf("ERROR: `u_b` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    if (w->stgs->do_super_scs == 0) {
        w->v = scs_calloc(l, sizeof (scs_float));
        if (w->v == SCS_NULL) {
            scs_printf("ERROR: `v` memory allocation failure\n");
            RETURN SCS_NULL;
        }
    }
    w->u_t = scs_malloc(l * sizeof (scs_float));
    if (w->u_t == SCS_NULL) {
        scs_printf("ERROR: `u_t` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->u_prev = scs_malloc(l * sizeof (scs_float));
    if (w->u_prev == SCS_NULL) {
        scs_printf("ERROR: `u_prev` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->h = scs_malloc((l - 1) * sizeof (scs_float));
    if (w->h == SCS_NULL) {
        scs_printf("ERROR: `h` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->g = scs_malloc((l - 1) * sizeof (scs_float));
    if (w->g == SCS_NULL) {
        scs_printf("ERROR: `g` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->pr = scs_malloc(d->m * sizeof (scs_float));
    if (w->pr == SCS_NULL) {
        scs_printf("ERROR: `pr` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->dr = scs_malloc(d->n * sizeof (scs_float));
    if (w->dr == SCS_NULL) {
        scs_printf("ERROR: `dr` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->b = scs_malloc(d->m * sizeof (scs_float));
    if (w->b == SCS_NULL) {
        scs_printf("ERROR: `b` memory allocation failure\n");
        RETURN SCS_NULL;
    }
    w->c = scs_malloc(d->n * sizeof (scs_float));
    if (w->c == SCS_NULL) {
        scs_printf("ERROR: `c` memory allocation failure\n");
        RETURN SCS_NULL;
    }



    if (w->stgs->do_super_scs == 1) {
        /* -------------------------------------
         * Additional memory needs to be allocated 
         * in SuperSCS
         * ------------------------------------- */
        w->R = scs_calloc(l, sizeof (scs_float));
        if (w->R == SCS_NULL) {
            scs_printf("ERROR: `R` memory allocation failure\n");
            RETURN SCS_NULL;
        }
        w->R_prev = scs_calloc(l, sizeof (scs_float));
        if (w->R_prev == SCS_NULL) {
            scs_printf("ERROR: `R_prev` memory allocation failure\n");
            RETURN SCS_NULL;
        }
        w->dir = scs_malloc(l * sizeof (scs_float));
        if (w->dir == SCS_NULL) {
            scs_printf("ERROR: `dir` memory allocation failure\n");
            RETURN SCS_NULL;
        }
        w->dut = scs_malloc(l * sizeof (scs_float));
        if (w->dut == SCS_NULL) {
            scs_printf("ERROR: `dut` memory allocation failure\n");
            RETURN SCS_NULL;
        }
        w->s_b = scs_malloc(d->m * sizeof (scs_float));
        if (w->s_b == SCS_NULL) {
            scs_printf("ERROR: `s_b` memory allocation failure\n");
            RETURN SCS_NULL;
        }

        w->stepsize = 1.0;

        /* -------------------------------------
         * Restarted Broyden requires the allocation
         * of an (S,U)-cache.
         * ------------------------------------- */
        if ((w->stgs->direction == restarted_broyden || w->stgs->direction == restarted_broyden_v2)
                && w->stgs->memory > 0) {
            w->su_cache = initSUCache(w->stgs->memory, l);
            if (w->su_cache == SCS_NULL) {
                scs_printf("ERROR: `su_cache` memory allocation failure\n");
                RETURN SCS_NULL;
            }
        } else {
            w->su_cache = SCS_NULL;
        }

        /* -------------------------------------
         * Allocate memory for the full Broyden
         * method
         * ------------------------------------- */
        if (w->stgs->direction == full_broyden) {
            scs_int i;
            w->H = scs_malloc(l * l * sizeof (scs_float));
            if (w->H == SCS_NULL) {
                scs_printf("ERROR: `H` memory allocation failure\n");
                RETURN SCS_NULL;
            }
            /* H = I */
            for (i = 0; i < l; ++i) {
                w->H[i * (l + 1)] = 1.0;
            }
        } else {
            w->H = SCS_NULL;
        }

        w->Sk = scs_malloc(l * sizeof (scs_float));
        if (w->Sk == SCS_NULL) {
            scs_printf("ERROR: `Sk` memory allocation failure\n");
            RETURN SCS_NULL;
        }
        w->Yk = scs_malloc(l * sizeof (scs_float));
        if (w->Yk == SCS_NULL) {
            scs_printf("ERROR: `Yk` memory allocation failure\n");
            RETURN SCS_NULL;
        }

        if (w->stgs->ls > 0) {
            w->wu = scs_malloc(l * sizeof (scs_float));
            if (w->wu == SCS_NULL) {
                scs_printf("ERROR: `wu` memory allocation failure\n");
                RETURN SCS_NULL;
            }
            w->Rwu = scs_malloc(l * sizeof (scs_float));
            if (w->Rwu == SCS_NULL) {
                scs_printf("ERROR: `Rwu` memory allocation failure\n");
                RETURN SCS_NULL;
            }
            w->wu_t = scs_malloc(l * sizeof (scs_float));
            if (w->wu_t == SCS_NULL) {
                scs_printf("ERROR: `wu_t` memory allocation failure\n");
                RETURN SCS_NULL;
            }
            w->wu_b = scs_malloc(l * sizeof (scs_float));
            if (w->wu_b == SCS_NULL) {
                scs_printf("ERROR: `wu_b` memory allocation failure\n");
                RETURN SCS_NULL;
            }
        }
    } else {
        /* -------------------------------------
         * In SCS, the pointers which correspond to
         * SuperSCS are set to SCS_NULL and are 
         * inactive.
         * ------------------------------------- */
        w->R = SCS_NULL;
        w->R_prev = SCS_NULL;
        w->dir = SCS_NULL;
        w->dut = SCS_NULL;
        w->s_b = SCS_NULL;
        w->su_cache = SCS_NULL;
        w->Yk = SCS_NULL;
        w->Sk = SCS_NULL;
        w->wu = SCS_NULL;
        w->Rwu = SCS_NULL;
        w->wu_t = SCS_NULL;
        w->wu_b = SCS_NULL;
    }

    w->A = d->A;
    if (w->stgs->normalize) {
#ifdef COPYAMATRIX
        if (!copyAMatrix(&(w->A), d->A)) {
            scs_printf("ERROR: copy A matrix failed\n");
            RETURN SCS_NULL;
        }
#endif
        w->scal = scs_malloc(sizeof (Scaling));
        normalizeA(w->A, w->stgs, k, w->scal);
#if EXTRAVERBOSE > 0
        printArray(w->scal->D, d->m, "D");
        scs_printf("norm D = %4f\n", calcNorm(w->scal->D, d->m));
        printArray(w->scal->E, d->n, "E");
        scs_printf("norm E = %4f\n", calcNorm(w->scal->E, d->n));
#endif
    } else {
        w->scal = SCS_NULL;
    }
    if (!(w->coneWork = initCone(k))) {
        scs_printf("ERROR: initCone failure\n");
        RETURN SCS_NULL;
    }
    w->p = initPriv(w->A, w->stgs);
    if (!w->p) {
        scs_printf("ERROR: initPriv failure\n");
        RETURN SCS_NULL;
    }
    RETURN w;
}

static scs_int updateWork(const Data *d, Work *w, const Sol *sol) {
    DEBUG_FUNC
    /* before normalization */
    scs_int n = d->n;
    scs_int m = d->m;

    w->nm_b = calcNorm(d->b, m);
    w->nm_c = calcNorm(d->c, n);
    memcpy(w->b, d->b, d->m * sizeof (scs_float));
    memcpy(w->c, d->c, d->n * sizeof (scs_float));

#if EXTRAVERBOSE > 0
    printArray(w->b, m, "b");
    scs_printf("pre-normalized norm b = %4f\n", calcNorm(w->b, m));
    printArray(w->c, n, "c");
    scs_printf("pre-normalized norm c = %4f\n", calcNorm(w->c, n));
#endif
    if (w->stgs->normalize) {
        normalizeBC(w);
#if EXTRAVERBOSE > 0
        printArray(w->b, m, "bn");
        scs_printf("sc_b = %4f\n", w->sc_b);
        scs_printf("post-normalized norm b = %4f\n", calcNorm(w->b, m));
        printArray(w->c, n, "cn");
        scs_printf("sc_c = %4f\n", w->sc_c);
        scs_printf("post-normalized norm c = %4f\n", calcNorm(w->c, n));
#endif
    }
    if (w->stgs->warm_start) {
        warmStartVars(w, sol);
    } else {
        coldStartVars(w);
    }
    memcpy(w->h, w->c, n * sizeof (scs_float));
    memcpy(&(w->h[n]), w->b, m * sizeof (scs_float));
    memcpy(w->g, w->h, (n + m) * sizeof (scs_float));
    solveLinSys(w->A, w->stgs, w->p, w->g, SCS_NULL, -1);
    scaleArray(&(w->g[n]), -1, m);
    w->gTh = innerProd(w->h, w->g, n + m);
    RETURN 0;
}

scs_int scs_solve(Work *w, const Data *d, const Cone *k, Sol *sol, Info *info) {
    DEBUG_FUNC
    scs_int i;
    timer solveTimer;
    struct residuals r;
    if (d == SCS_NULL
            || k == SCS_NULL
            || sol == SCS_NULL
            || info == SCS_NULL
            || w == SCS_NULL
            || d->b == SCS_NULL
            || d->c == SCS_NULL) {
        scs_printf("ERROR: SCS_NULL input\n");
        RETURN SCS_FAILED;
    }
    /* initialize ctrl-c support */
    startInterruptListener();
    tic(&solveTimer);
    info->statusVal = SCS_UNFINISHED; /* not yet converged */
    r.lastIter = -1;
    updateWork(d, w, sol);

    if (w->stgs->verbose)
        printHeader(w, k);
    /* scs: */
    for (i = 0; i < w->stgs->max_iters; ++i) {
        memcpy(w->u_prev, w->u, w->l * sizeof (scs_float));

        if (projectLinSys(w, i) < 0) {
            RETURN failure(w, w->m, w->n, sol, info, SCS_FAILED,
                    "error in projectLinSys", "Failure");
        }
        if (projectCones(w, k, i) < 0) {
            RETURN failure(w, w->m, w->n, sol, info, SCS_FAILED,
                    "error in projectCones", "Failure");
        }

        updateDualVars(w);

        if (isInterrupted()) {
            RETURN failure(w, w->m, w->n, sol, info, SCS_SIGINT, "Interrupted",
                    "Interrupted");
        }
        if (i % CONVERGED_INTERVAL == 0) {
            calcResiduals(w, &r, i);
            if ((info->statusVal = hasConverged(w, &r, i)) != 0) {
                break;
            }
        }

        if (w->stgs->verbose && i % PRINT_INTERVAL == 0) {
            calcResiduals(w, &r, i);
            printSummary(w, i, &r, &solveTimer);
        }
    }
    if (w->stgs->verbose) {
        calcResiduals(w, &r, i);
        printSummary(w, i, &r, &solveTimer);
    }
    /* populate solution vectors (unnormalized) and info */
    getSolution(w, sol, info, &r, i);
    info->solveTime = tocq(&solveTimer);

    if (w->stgs->verbose)
        printFooter(d, k, sol, w, info); /* LCOV_EXCL_LINE */
    endInterruptListener();
    RETURN info->statusVal;
}

static void compute_sb_kapb(scs_float *u, scs_float *u_b, scs_float *u_t, Work *w) {
    scs_int j;
    for (j = 0; j < w->m; ++j) {
        w->s_b[j] = u_b[w->n + j] - 2.0 * u_t[w->n + j] + u[w->n + j];
    }
    w->kap_b = u_b[w->l - 1] - 2.0 * u_t[w->l - 1] + u[w->l - 1];
}

scs_int superscs_solve(Work *work, const Data *data, const Cone *cone, Sol *sol, Info *info) {
    DEBUG_FUNC
    scs_int i; /* i indexes the (outer) iterations */
    scs_int j; /* j indexes the line search iterations */
    scs_int j1; /* j1 indexes other auxiliary iterations (e.g., algebraic operations) */
    scs_int how = 0; /* -1:unsuccessful backtracking, 0:K0, 1:K1, 2:K2 */
    scs_float eta;
    scs_float r_safe;
    scs_float nrmRw_con; /* norm of FP res at line-search */
    scs_float nrmR_con_old; /* keeps previous FP res */
    scs_float slack; /* for K2 */
    scs_float rhs; /* for K2 */
    scs_float stepsize2; /* for K2 */
    scs_float sqrt_rhox = SQRTF(work->stgs->rho_x);
    scs_float q = work->stgs->sse;
    timer solveTimer;
    struct residuals r;

    if (work->stgs->do_record_progress) {
        /*TODO store the norm of FPR */
        const scs_int max_history_alloc = data->stgs->max_iters / CONVERGED_INTERVAL;
        info->progress_relgap = malloc(sizeof (scs_float) * max_history_alloc);
        info->progress_respri = malloc(sizeof (scs_float) * max_history_alloc);
        info->progress_resdual = malloc(sizeof (scs_float) * max_history_alloc);
        info->progress_pcost = malloc(sizeof (scs_float) * max_history_alloc);
        info->progress_dcost = malloc(sizeof (scs_float) * max_history_alloc);
        info->progress_iter = malloc(sizeof (scs_int) * max_history_alloc);
        info->progress_norm_fpr = malloc(sizeof (scs_int) * max_history_alloc);
    } else {
        info->progress_relgap = SCS_NULL;
        info->progress_respri = SCS_NULL;
        info->progress_resdual = SCS_NULL;
        info->progress_pcost = SCS_NULL;
        info->progress_dcost = SCS_NULL;
        info->progress_iter = SCS_NULL;
        info->progress_norm_fpr = SCS_NULL;
    }

    if (data == SCS_NULL
            || cone == SCS_NULL
            || sol == SCS_NULL
            || info == SCS_NULL
            || work == SCS_NULL
            || data->b == SCS_NULL
            || data->c == SCS_NULL) {
        scs_printf("ERROR: SCS_NULL input\n");
        RETURN SCS_FAILED;
    }

    /* initialize ctrl-c support */
    startInterruptListener();
    tic(&solveTimer);
    info->statusVal = SCS_UNFINISHED; /* not yet converged */
    r.lastIter = -1;
    updateWork(data, work, sol);

    if (work->stgs->verbose > 0)
        printHeader(work, cone);

    /* Initialize: */
    i = 0; /* Needed for the next two functions */
    if (projectLinSysv2(work->u_t, work->u, work, i) < 0) { /* u_t = (I+Q)^{-1} u*/
        RETURN failure(work, work->m, work->n, sol, info, SCS_FAILED,
                "error in projectLinSysv2", "Failure");
    }
    if (projectConesv2(work->u_b, work->u_t, work->u, work, cone, i) < 0) { /* u_bar = proj_C(2u_t - u) */
        RETURN failure(work, work->m, work->n, sol, info, SCS_FAILED,
                "error in projectConesv2", "Failure");
    }
    compute_sb_kapb(work->u, work->u_b, work->u_t, work); /* compute s_b and kappa_b */
    calcFPRes(work->R, work->u_t, work->u_b, work->l); /* compute Ru */
    scaleArray(work->R, sqrt_rhox, work->n); /* scale R_x with sqrt_rhox */
    eta = calcNorm(work->R, work->l); /* initialize eta = |Ru^0| (norm of scaled R) */
    scaleArray(work->u, sqrt_rhox, work->n); /* u is now scaled */
    r_safe = eta;

    work->nrmR_con = eta;

    /***** HENCEFORTH, R IS SCALED! *****/

    /* MAIN SUPER SCS LOOP */
    for (i = 0; i < work->stgs->max_iters; ++i) {

        if (isInterrupted()) {
            RETURN failure(work, work->m, work->n, sol, info, SCS_SIGINT, "Interrupted",
                    "Interrupted");
        }

        /* Convergence checks */
        if (i % CONVERGED_INTERVAL == 0) {
            calcResidualsSuperscs(work, &r, i);
            if (work->stgs->do_record_progress) {
                scs_int idx_progress = i / CONVERGED_INTERVAL;
                info->progress_iter[idx_progress] = i;
                info->progress_relgap[idx_progress] = r.relGap;
                info->progress_respri[idx_progress] = r.resPri;
                info->progress_resdual[idx_progress] = r.resDual;
                info->progress_pcost[idx_progress] = r.cTx_by_tau / r.tau;
                info->progress_dcost[idx_progress] = -r.bTy_by_tau / r.tau;
            }
            if ((info->statusVal = hasConverged(work, &r, i))) {
                break;
            }
        }

        /* Prints results every PRINT_INTERVAL iterations */
        if (work->stgs->verbose && i % PRINT_INTERVAL == 0) {
            calcResidualsSuperscs(work, &r, i);
            printSummary(work, i, &r, &solveTimer);
        }

        if (work->stgs->ls > 0 || work->stgs->k0 == 1) {
            work->stgs->sse *= q; /*sse = q^i */
            if (i == 0) {
                /* -------------------------------------------
                 * At i=0, the direction is defined using the 
                 * FPR: dir^0 = -R 
                 * -------------------------------------------- */
                setAsScaledArray(work->dir, work->R, -1, work->l);
            } else {
                if (how == 0 || work->stgs->ls == 0) {
                    for (j1 = 0; j1 < work->l; ++j1) {
                        work->Sk[j1] = work->u[j1] - work->u_prev[j1];
                        work->Yk[j1] = work->R[j1] - work->R_prev[j1];
                    }
                } else {
                    for (j1 = 0; j1 < work->l; ++j1) {
                        work->Sk[j1] = work->wu[j1] - work->u_prev[j1];
                        work->Yk[j1] = work->Rwu[j1] - work->R_prev[j1];
                    }
                }
                /* compute direction */
                computeDirection(work, i);
            }
            /* -------------------------------------------
             * Scale the x-part of dir using sqrt_rhox
             * -------------------------------------------- */
            scaleArray(work->dir, 1 / sqrt_rhox, work->n);
        }

        memcpy(work->u_prev, work->u, work->l * sizeof (scs_float)); /* u_prev = u */
        memcpy(work->R_prev, work->R, work->l * sizeof (scs_float)); /* R_prev = R */
        how = -1; /* no backtracking (yet) */
        nrmR_con_old = work->nrmR_con;

        if (i >= work->stgs->warm_start) {
            if (work->stgs->k0 == 1 && work->nrmR_con <= work->stgs->c_bl * eta) {
                addArray(work->u, work->dir, work->l); /* u += dir */
                how = 0;
                eta = work->nrmR_con;
                work->stepsize = 1.0;
            } else if (work->stgs->ls > 0) {
                projectLinSysv2(work->dut, work->dir, work, i);
                work->stepsize = 2.0;

                /* Line search */
                for (j = 0; j < work->stgs->ls; ++j) {
                    work->stepsize *= work->stgs->beta;
                    for (j1 = 0; j1 < work->l; ++j1) {
                        work->wu[j1] = work->u[j1] + work->stepsize * work->dir[j1]; /* wu = u + step * dir */
                        work->wu_t[j1] = work->u_t[j1] + work->stepsize * work->dut[j1]; /* wut = u_t + step * dut */
                    }

                    projectConesv2(work->wu_b, work->wu_t, work->wu, work, cone, i);
                    calcFPRes(work->Rwu, work->wu_t, work->wu_b, work->l); /* calculate FPR on scaled vectors */

                    nrmRw_con = calcNorm(work->Rwu, work->l);
                    /* K1 */
                    if (work->stgs->k1
                            && nrmRw_con <= work->stgs->c1 * nrmR_con_old
                            && work->nrmR_con <= r_safe) { /* a bit different than matlab */
                        memcpy(work->u, work->wu, work->l * sizeof (scs_float));
                        memcpy(work->u_t, work->wu_t, work->l * sizeof (scs_float));
                        memcpy(work->u_b, work->wu_b, work->l * sizeof (scs_float));
                        memcpy(work->R, work->Rwu, work->l * sizeof (scs_float));
                        compute_sb_kapb(work->wu, work->wu_b, work->wu_t, work);
                        work->nrmR_con = nrmRw_con;
                        r_safe = work->nrmR_con + work->stgs->sse; /* The power already computed at the beginning of the main loop */
                        how = 1;
                        break;
                    }

                    /* K2 */
                    if (work->stgs->k2) {
                        slack = nrmRw_con * nrmRw_con - work->stepsize * innerProd(work->dir, work->Rwu, work->l);
                        rhs = work->stgs->sigma * work->nrmR_con * nrmRw_con;
                        /* printf("%2.12f, %2.12f \n", slack, rhs); */
                        if (slack >= rhs) {
                            stepsize2 = (work->stgs->alpha * (slack / (nrmRw_con * nrmRw_con)));
                            addScaledArray(work->u, work->Rwu, work->l, -stepsize2);
                            how = 2;
                            break; /* exits the line search loop */
                        }
                    }

                } /* end of line-search */
            }

        }
        if (how == -1) { /* means that R didn't change */
            /* x -= alpha*sqrt(rho)*Rx */
            addScaledArray(work->u, work->R, work->n, -work->stgs->alpha * sqrt_rhox);
            addScaledArray(work->u + work->n, work->R + work->n, work->m + 1, -work->stgs->alpha);
        }
        if (how != 1) { /* exited with K1 */
            projectLinSysv2(work->u_t, work->u, work, i);
            projectConesv2(work->u_b, work->u_t, work->u, work, cone, i); /* u_bar = proj_C(2u_t - u) */
            compute_sb_kapb(work->u, work->u_b, work->u_t, work);
            calcFPRes(work->R, work->u_t, work->u_b, work->l);
            scaleArray(work->R, sqrt_rhox, work->n);
            work->nrmR_con = calcNorm(work->R, work->l);
            /* printf("%3.12f\n", work->nrmR_con); */
        }
    } /* main for loop */

    /* prints summary of last iteration */
    if (work->stgs->verbose) {
        calcResidualsSuperscs(work, &r, i);
        printSummary(work, i, &r, &solveTimer);
    }

    /* populate solution vectors (unnormalized) and info */
    /* update u, denormalize, etc */
    getSolution(work, sol, info, &r, i);
    info->iter = i;
    info->solveTime = tocq(&solveTimer);
    info->history_length = i / CONVERGED_INTERVAL;

    if (work->stgs->verbose)
        printFooter(data, cone, sol, work, info); /* LCOV_EXCL_LINE */
    endInterruptListener();

    RETURN info->statusVal;
}

void scs_finish(Work * w) {
    DEBUG_FUNC
    if (w) {
        finishCone(w->coneWork);
        if (w->stgs && w->stgs->normalize) {
#ifndef COPYAMATRIX
            unNormalizeA(w->A, w->stgs, w->scal);
#else
            freeAMatrix(w->A);
#endif
        }

        if (w->p)
            freePriv(w->p);
        freeWork(w);
    }
#if EXTRAVERBOSE > 0
    scs_printf("exit finish\n");
#endif
    RETURN;
}

Work * scs_init(const Data *d, const Cone *k, Info * info) {
    DEBUG_FUNC
#if EXTRAVERBOSE > 1
            tic(&globalTimer);
#endif
    Work *w;
    timer initTimer;
    startInterruptListener();
    if (d == SCS_NULL
            || k == SCS_NULL
            || info == SCS_NULL) {
        /* LCOV_EXCL_START */
        scs_printf("ERROR: Missing Data, Cone or Info input\n");
        RETURN SCS_NULL;
        /* LCOV_EXCL_STOP */
    }
#if EXTRAVERBOSE > 0
    printData(d);
    printConeData(k);
#endif
#ifndef NOVALIDATE
    if (validate(d, k) < 0) {
        scs_printf("ERROR: Validation returned failure\n");
        RETURN SCS_NULL;
    }
#endif
    tic(&initTimer);
    w = initWork(d, k);
    /* strtoc("init", &initTimer); */
    info->setupTime = tocq(&initTimer);
    if (d->stgs->verbose) {

        scs_printf("Setup time: %1.2es\n", info->setupTime / 1e3);
    }
    endInterruptListener();
    RETURN w;
}

/* this just calls scs_init, scs_solve, and scs_finish */
scs_int scs(const Data *d, const Cone *k, Sol *sol, Info * info) {

    DEBUG_FUNC
    scs_int status;
    Work *w = scs_init(d, k, info);

    if (d->stgs->verbose >= 2) {
        /* LCOV_EXCL_START */
        scs_printf("Settings:\n"
                "alpha        : %g\n"
                "beta         : %g\n"
                "c1           : %g\n"
                "c_bl         : %g\n"
                "cg_rate      : %g\n"
                "dir          : %d\n"
                "do_super_scs : %d\n"
                "eps          : %g\n"
                "(k0, k1, k2) : (%d, %d, %d)\n"
                "ls           : %d\n"
                "max_iters    : %d\n"
                "memory       : %d\n"
                "normalize    : %d\n"
                "rho_x        : %g\n"
                "scale        : %g\n"
                "sigma        : %g\n"
                "sse          : %g\n"
                "thetabar     : %g\n"
                "warm_start   : %d\n",
                d->stgs->alpha,
                d->stgs->beta,
                d->stgs->c1,
                d->stgs->c_bl,
                d->stgs->cg_rate,
                d->stgs->direction,
                d->stgs->do_super_scs,
                d->stgs->eps,
                d->stgs->k0,
                d->stgs->k1,
                d->stgs->k2,
                d->stgs->ls,
                d->stgs->max_iters,
                d->stgs->memory,
                d->stgs->normalize,
                d->stgs->rho_x,
                d->stgs->scale,
                d->stgs->sigma,
                d->stgs->sse,
                d->stgs->thetabar,
                d->stgs->warm_start);
        /* LCOV_EXCL_STOP */
    }

#if EXTRAVERBOSE > 0
    scs_printf("size of scs_int = %lu, size of scs_float = %lu\n",
            sizeof (scs_int), sizeof (scs_float));
#endif
    if (w) {
        if (w->stgs->do_super_scs) {
            /* solve with SuperSCS*/
            superscs_solve(w, d, k, sol, info);
        } else {
            /* solve with SCS */
            scs_solve(w, d, k, sol, info);
        }
        status = info->statusVal;
    } else {
        status = failure(SCS_NULL, d ? d->m : -1, d ? d->n : -1, sol, info,
                SCS_FAILED, "could not initialize work", "Failure");
    }
    scs_finish(w);
    RETURN status;
}

Sol * initSol() {
    Sol *sol = scs_calloc(1, sizeof (* sol));
    if (sol == SCS_NULL) {
        /* LCOV_EXCL_START */
        scs_printf("ERROR: allocating sol failure\n");
        RETURN SCS_NULL;
        /* LCOV_EXCL_STOP */
    }
    sol->s = SCS_NULL;
    sol->x = SCS_NULL;
    sol->y = SCS_NULL;
    RETURN sol;
}

Info * initInfo() {
    Info * info = scs_calloc(1, sizeof (*info));
    if (info == SCS_NULL) {
        /* LCOV_EXCL_START */
        scs_printf("ERROR: allocating info failure\n");
        RETURN SCS_NULL;
        /* LCOV_EXCL_STOP */
    }
    info->pobj = NAN;
    info->dobj = NAN;
    info->iter = -1;
    info->relGap = NAN;
    info->resDual = NAN;
    info->resInfeas = NAN;
    info->resPri = NAN;
    info->resUnbdd = NAN;
    info->setupTime = NAN;
    info->solveTime = NAN;
    info->statusVal = SCS_INDETERMINATE;
    info->progress_iter = SCS_NULL;
    info->progress_dcost = SCS_NULL;
    info->progress_pcost = SCS_NULL;
    info->progress_relgap = SCS_NULL;
    info->progress_resdual = SCS_NULL;
    info->progress_respri = SCS_NULL;
    RETURN info;
}

Data * initData() {
    Data * data = malloc(sizeof (*data));

    if (data == SCS_NULL) {
        /* LCOV_EXCL_START */
        scs_printf("ERROR: allocating data failure\n");
        RETURN SCS_NULL;
        /* LCOV_EXCL_STOP */
    }

    data->A = SCS_NULL;
    data->b = SCS_NULL;
    data->c = SCS_NULL;
    data->m = 0;
    data->n = 0;

    data->stgs = scs_malloc(sizeof (Settings));
    setDefaultSettings(data);

    RETURN data;
}
