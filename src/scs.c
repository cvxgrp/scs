#include "scs.h"
#include "normalize.h"
#include "directions.h"

#ifndef EXTRAVERBOSE
/* if verbose print summary output every this num iterations */
#define PRINT_INTERVAL 100
/* check for convergence every this num iterations */
#define CONVERGED_INTERVAL 20
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
    " pri obj ", " dua obj ", " kap/tau ", " time (s)",
};
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 8;
static const scs_int LINE_LEN = 76;

static scs_int scs_isnan(scs_float x) {
    DEBUG_FUNC
    RETURN(x == NAN || x != x);
}

static SUCache * initSUCache(scs_int memory, scs_int l) {
    SUCache * cache = scs_calloc(1, sizeof (SUCache));
    if (!cache) {
        scs_printf("ERROR: allocating YSCache failure\n");
        RETURN SCS_NULL;
    }

    /* we allocate one extra memory position because it's needed */
    cache->S = scs_calloc((1 + memory) * l, sizeof (scs_float)); /* S: l-by-mem */
    cache->U = scs_calloc((1 + memory) * l, sizeof (scs_float)); /* U: l-by-mem */


    /* the cache must know its memory length */
    cache->mem = memory;

    /* initial active memory is 0 */
    resetSUCache(cache);
    return cache;
}

static void freeYSCache(SUCache * cache) {
    if (!cache) {
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
    if (!w)
        RETURN;
    if (w->u)
        scs_free(w->u);
    if (w->v)
        scs_free(w->v);
    if (w->u_t)
        scs_free(w->u_t);
    if (w->u_prev)
        scs_free(w->u_prev);
    if (w->h)
        scs_free(w->h);
    if (w->g)
        scs_free(w->g);
    if (w->b)
        scs_free(w->b);
    if (w->c)
        scs_free(w->c);
    if (w->pr)
        scs_free(w->pr);
    if (w->dr)
        scs_free(w->dr);
    if (w->scal) {
        if (w->scal->D)
            scs_free(w->scal->D);
        if (w->scal->E)
            scs_free(w->scal->E);
        scs_free(w->scal);
    }
    if (w->u_b)
        scs_free(w->u_b);
    if (w->R)
        scs_free(w->R);
    if (w->R_prev)
        scs_free(w->R_prev);
    if (w->dir)
        scs_free(w->dir);
    if (w->dut)
        scs_free(w->dut);
    if (w->Sk)
        scs_free(w->Sk);
    if (w->Yk)
        scs_free(w->Yk);
    if (w->wu)
        scs_free(w->wu);
    if (w->wu_t)
        scs_free(w->wu_t);
    if (w->wu_b)
        scs_free(w->wu_b);
    if (w->Rwu)
        scs_free(w->Rwu);
    if (w->su_cache)
        freeYSCache(w->su_cache);

    scs_free(w);
    RETURN;
}

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
            if (!sol->x)
                sol->x = scs_malloc(sizeof (scs_float) * n);
            scaleArray(sol->x, NAN, n);
        }
        if (m > 0) {
            if (!sol->y)
                sol->y = scs_malloc(sizeof (scs_float) * m);
            scaleArray(sol->y, NAN, m);
            if (!sol->s)
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
        if (!w->stgs->do_super_scs && scs_isnan(w->v[i])) {
            w->v[i] = 0;
        }
    }
#endif
    if (w->stgs->normalize) {
        normalizeWarmStart(w);
    }
    RETURN;
}

static void warmStartVarsv2(Work *w, const Sol *sol) {
    DEBUG_FUNC
    scs_int i, n = w->n, m = w->m;

    memcpy(w->u, sol->x, n * sizeof (scs_float));
    memcpy(&(w->u[n]), sol->y, m * sizeof (scs_float));
    w->u[n + m] = 1.0;
#ifndef NOVALIDATE
    for (i = 0; i < n + m + 1; ++i) {
        if (scs_isnan(w->u[i]))
            w->u[i] = 0;
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
    RETURN SQRTF(pres); /* norm(Ax + s - b * tau) */
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
    RETURN SQRTF(dres); /* norm(A'y + c * tau) */
}

/* calculates un-normalized quantities */
static void calcResiduals(Work *w, struct residuals *r, scs_int iter) {
    DEBUG_FUNC
    scs_float *x = w->u, *y = &(w->u[w->n]), *s = &(w->v[w->n]);
    scs_float nmpr_tau, nmdr_tau, nmAxs_tau, nmATy_tau, cTx, bTy;
    scs_int n = w->n, m = w->m;

    /* checks if the residuals are unchanged by checking iteration */
    if (r->lastIter == iter) {
        RETURN;
    }
    r->lastIter = iter;

    r->tau = ABS(w->u[n + m]);
    r->kap = ABS(w->v[n + m]) /
            (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);

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
    DEBUG_FUNC
            const scs_int l = w->l;
    scs_int status;

    /* u_t = u (is already provided scaled      */
    memcpy(u_t, u, l * sizeof (scs_float));

    /* ut(1:l-1) = ut(1:l-1) - ut(end) * h      */
    addScaledArray(u_t, w->h, l - 1, -u_t[l - 1]);

    /* ut -= scalar * h                         */
    addScaledArray(u_t, w->h, l - 1,
            -innerProd(u_t, w->g, l - 1) / (w->gTh + 1));

    /* ut(n+1:end-1) = -ut(n+1:end-1);           */
    scaleArray(u_t + w->n, -1, w->m);

    /* call `solveLinSys` to update ut(1:n+m)   */
    status = solveLinSys(w->A, w->stgs, w->p, u_t, u, iter);

    /* ut(end) = (ut(end) + h'*ut(1:l-1))       */
    u_t[l - 1] += innerProd(u_t, w->h, l - 1);

    RETURN status;
}

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
    if (u_b[l - 1] < 0.0)
        u_b[l - 1] = 0.0;

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
    memcpy(sol->y, &(w->u[w->n]), w->m * sizeof (scs_float));
    RETURN;
}

static void sets(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->s) {
        sol->s = scs_malloc(sizeof (scs_float) * w->m);
    }
    memcpy(sol->s, &(w->v[w->n]), w->m * sizeof (scs_float));
    RETURN;
}

static void setx(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->x)
        sol->x = scs_malloc(sizeof (scs_float) * w->n);
    memcpy(sol->x, w->u, w->n * sizeof (scs_float));
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
static void getSolution(Work *w, Sol *sol, Info *info, struct residuals *r,
        scs_int iter) {
    DEBUG_FUNC
    scs_int l = w->n + w->m + 1;
    calcResiduals(w, r, iter);
    setx(w, sol);
    sety(w, sol);
    sets(w, sol);
    if (info->statusVal == SCS_UNFINISHED) {
        /* not yet converged, take best guess */
        if (r->tau > INDETERMINATE_TOL && r->tau > r->kap) {
            info->statusVal = solved(w, sol, info, r->tau);
        } else if (calcNorm(w->u, l) <
                INDETERMINATE_TOL * SQRTF((scs_float) l)) {
            info->statusVal = indeterminate(w, sol, info);
        } else if (r->bTy_by_tau < r->cTx_by_tau) {
            info->statusVal = infeasible(w, sol, info, r->bTy_by_tau);
        } else {
            info->statusVal = unbounded(w, sol, info, r->cTx_by_tau);
        }
    } else if (isSolvedStatus(info->statusVal)) {
        info->statusVal = solved(w, sol, info, r->tau);
    } else if (isInfeasibleStatus(info->statusVal)) {
        info->statusVal = infeasible(w, sol, info, r->bTy_by_tau);
    } else {
        info->statusVal = unbounded(w, sol, info, r->cTx_by_tau);
    }
    if (w->stgs->normalize) {
        unNormalizeSol(w, sol);
    }
    getInfo(w, sol, info, r, iter);
    RETURN;
}

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
    scs_printf("%*.2e ", (int) HSPACE, tocq(solveTimer) / 1e3);
    scs_printf("\n");

#if EXTRAVERBOSE > 0
    scs_printf("Norm u = %4f, ", calcNorm(w->u, w->n + w->m + 1));
    scs_printf("Norm u_t = %4f, ", calcNorm(w->u_t, w->n + w->m + 1));
    scs_printf("Norm v = %4f, ", calcNorm(w->v, w->n + w->m + 1));
    scs_printf("tau = %4f, ", w->u[w->n + w->m]);
    scs_printf("kappa = %4f, ", w->v[w->n + w->m]);
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
    for (i = 0; i < HEADER_LEN - 1; ++i) {
        scs_printf("%s|", HEADER[i]);
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
        scs_printf("m and n must both be greater than 0; m = %li, n = %li\n",
                (long) d->m, (long) d->n);
        RETURN - 1;
    }
    if (d->m < d->n) {
        scs_printf("WARN: m less than n, problem likely degenerate\n");
        /* RETURN -1; */
    }
    if (validateLinSys(d->A) < 0) {
        scs_printf("invalid linear system input data\n");
        RETURN - 1;
    }
    if (validateCones(d, k) < 0) {
        scs_printf("cone validation error\n");
        RETURN - 1;
    }
    if (stgs->max_iters <= 0) {
        scs_printf("max_iters must be positive (max_iters=%d)\n", stgs->max_iters);
        RETURN - 1;
    }
    if (stgs->eps <= 0) {
        scs_printf("eps tolerance must be positive (eps=%g)\n", stgs->eps);
        RETURN - 1;
    }
    if (stgs->alpha <= 0 || stgs->alpha >= 2) {
        scs_printf("alpha must be in (0,2) (alpha=%g)\n", stgs->alpha);
        RETURN - 1;
    }
    if (stgs->rho_x <= 0) {
        scs_printf("rhoX must be positive (1e-3 works well) (rho_x=%g).\n", stgs->rho_x);
        RETURN - 1;
    }
    if (stgs->scale <= 0) {
        scs_printf("Parameter `scale` must be positive (1 works well).\n");
        RETURN - 1;
    }
    /* validate settings related to SuperSCS */
    if (stgs->do_super_scs == 1) {
        if (stgs->thetabar < 0 || stgs->thetabar > 1) {
            scs_printf("Parameters `thetabar` must be a scalar between 0 and 1 (thetabar=%g)\n", stgs->thetabar);
            RETURN - 1;
        }
        if (stgs->memory <= 1) {
            scs_printf("Quasi-Newton memory length (mem=%d) is too low; choose an integer at least equal to 2.\n", stgs->memory);
            RETURN - 1;
        }
        if (stgs->beta >= 1 || stgs->beta <= 0) {
            scs_printf("Stepsize reduction factor (beta=%g) out of bounds.\n", stgs->beta);
            RETURN - 1;
        }
        if (stgs->ls < 0) {
            scs_printf("Illegal maximum number of line search iterations (ls=%d).\n", stgs->ls);
            RETURN - 1;
        }
        if (stgs->sigma < 0) {
            scs_printf("Parameter sigma of the line search (sigma=%g) cannot be negative.\n", stgs->sigma);
            RETURN - 1;
        }
        if (stgs->c_bl < 0 || stgs->c1 >= 1) {
            scs_printf("Parameter (c_0=%g) for blind updates out of bounds.\n", stgs->c_bl);
            RETURN - 1;
        }
        if (stgs->c1 < 0 || stgs->c1 >= 1) {
            scs_printf("Parameter (c1=%g) for step K1 out of bounds.\n", stgs->c1);
            RETURN - 1;
        }
        if (stgs->sse < 0 || stgs->sse >= 1) {
            scs_printf("Parameter (sse=%g) for step K1 out of bounds.\n", stgs->sse);
            RETURN - 1;
        }
    }
    RETURN 0;
}

static Work *initWork(const Data *d, const Cone *k) {
    DEBUG_FUNC
    Work *w = scs_calloc(1, sizeof (Work));
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

    /* allocate workspace: */
    w->u = scs_malloc(l * sizeof (scs_float));
    w->v = scs_malloc(l * sizeof (scs_float));
    w->u_t = scs_malloc(l * sizeof (scs_float));
    w->u_prev = scs_malloc(l * sizeof (scs_float));
    w->h = scs_malloc((l - 1) * sizeof (scs_float));
    w->g = scs_malloc((l - 1) * sizeof (scs_float));
    w->pr = scs_malloc(d->m * sizeof (scs_float));
    w->dr = scs_malloc(d->n * sizeof (scs_float));
    w->b = scs_malloc(d->m * sizeof (scs_float));
    w->c = scs_malloc(d->n * sizeof (scs_float));

    /* added for superscs*/
    w->u_t = scs_malloc(l * sizeof (scs_float));
    w->R = scs_calloc(l, sizeof (scs_float));
    w->R_prev = scs_calloc(l, sizeof (scs_float));
    w->dir = scs_calloc(l, sizeof (scs_float));
    w->dut = scs_calloc(l, sizeof (scs_float));

    w->stepsize = 1.0;

    /* make cache */
    if (w->stgs->memory > 0) {
        w->su_cache = initSUCache(w->stgs->memory, l);
    } else {
        w->su_cache = SCS_NULL;
    }

    w->Sk = scs_calloc(l, sizeof (scs_float));
    w->Yk = scs_calloc(l, sizeof (scs_float));

    if (w->stgs->ls > 0) {
        w->wu = scs_calloc(l, sizeof (scs_float));
        w->Rwu = scs_calloc(l, sizeof (scs_float));
        w->wu_t = scs_calloc(l, sizeof (scs_float));
        w->wu_b = scs_calloc(l, sizeof (scs_float));
    }

    if (!w->u || !w->v || !w->u_t || !w->u_prev || !w->h || !w->g || !w->pr ||
            !w->dr || !w->b || !w->c || !w->u_b || !w->R || !w->R_prev ||
            !w->dir || !w->dut || !w->Sk || !w->Yk
            || (w->stgs->ls > 0 && (!w->wu || !w->Rwu || !w->wu_t || !w->wu_b))) {
        scs_printf("ERROR: work memory allocation failure\n");
        RETURN SCS_NULL;
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
    if (!d || !k || !sol || !info || !w || !d->b || !d->c) {
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
        printFooter(d, k, sol, w, info);
    endInterruptListener();
    RETURN info->statusVal;
}

scs_int superscs_solve(Work *w, const Data *d, const Cone *k, Sol *sol, Info *info) {
    DEBUG_FUNC
    scs_int i;
    scs_int j;
    scs_int how = 0; /* -1:unsuccessful backtracking, 0:K0, 1:K1, 2:K2 */
    scs_float eta;
    scs_float r_safe;
    scs_float nrmRw_con; /* norm of FP res at line-search */
    scs_float nrmR_con_old; /* keeps previous FP res */
    scs_float slack; /* for K2 */
    scs_float rhs; /* for K2 */
    scs_float stepsize2; /* for K2 */

    timer solveTimer;
    struct residuals r;
    if (!d || !k || !sol || !info || !w || !d->b || !d->c) {
        scs_printf("ERROR: SCS_NULL input\n");
        RETURN SCS_FAILED;
    }

    /* initialize ctrl-c support */
    startInterruptListener();
    tic(&solveTimer);
    info->statusVal = SCS_UNFINISHED; /* not yet converged */
    r.lastIter = -1;
    /*  updateWorkv2(d, w, sol); */

    if (w->stgs->verbose)
        printHeader(w, k);

    /* Initialize: */
    i = 0; /* Needed for the next two functions */
    projectLinSysv2(w->u_t, w->u, w, i); /* u_t = (I+Q)^{-1} u*/
    projectConesv2(w->u_b, w->u_t, w->u, w, k, i); /* u_bar = proj_C(2u_t - u) */

    calcFPRes(w->R, w->u_t, w->u_b, w->l); /* computes Ru */
    scaleArray(w->R, sqrt(w->stgs->rho_x), w->n);
    eta = calcNorm(w->R, w->l); /* initialize eta = |Ru^0| (norm of scaled R) */
    scaleArray(w->u, sqrt(w->stgs->rho_x), w->n); /* u is now scaled */
    r_safe = eta;

    /***** HENCEFORTH, R and u ARE SCALED! *****/

    /* MAIN SUPER SCS LOOP */
    for (i = 0; i < w->stgs->max_iters; ++i) {

        w->nrmR_con = (i == 0) ? eta : calcNorm(w->R, w->l);

        if (w->stgs->ls > 0 || w->stgs->k0 == 1) {
            if (i == 0) {
                memcpy(w->dir, w->R, w->l * sizeof(scs_float));
                scaleArray(w->dir, -1, w->l); /* dir^0 = -R */
            } else {
                w->stgs->sse *= w->stgs->sse;
                /*  if (w->how[i-1] == 0 || w->stgs->ls == 0) { */
                if (how == 0 || w->stgs->ls == 0) {
                    memcpy(w->Sk, w->u, w->l * sizeof(scs_float));
                    addScaledArray(w->Sk, w->u_prev, w->l, -1);
                    memcpy(w->Yk, w->R, w->l * sizeof(scs_float));
                    addScaledArray(w->Yk, w->R_prev, w->l, -1);
                } else {
                    memcpy(w->Sk, w->wu, w->l * sizeof(scs_float));
                    addScaledArray(w->Sk, w->u_prev, w->l, -1);
                    memcpy(w->Yk, w->Rwu, w->l * sizeof(scs_float));
                    addScaledArray(w->Yk, w->R_prev, w->l, -1);
                }
                /* compute direction */
                computeLSBroyden(w);
                scaleArray(w->dir, 1 / sqrt(w->stgs->rho_x), w->n);
            }
        }
        memcpy(w->u_prev, w->u, w->l * sizeof (scs_float));
        how = -1; /* no backtracking (yet) */
        nrmR_con_old = w->nrmR_con;
        if (i >= w->stgs->warm_start) {
            if (w->stgs->k0 == 1 && w->nrmR_con <= w->stgs->c_bl * eta) {
                addScaledArray(w->u, w->dir, w->l, 1.0);
                how = 0;
                eta = w->nrmR_con;
                r_safe = INFINITY; // TODO: chk if it should be inf.
                w->stepsize = 1.0;
            } else if (w->stgs->ls > 0) {
                projectLinSysv2(w->dut, w->dir, w, i);
                w->stepsize = 2.0;

                /* Line - search */
                for (j = 0; j < w->stgs->ls; ++j) {
                    w->stepsize *= w->stgs->beta;
                    addScaledArray(w->wu, w->dir, w->l, w->stepsize);
                    addScaledArray(w->wu_t, w->dut, w->l, w->stepsize);
                    projectConesv2(w->wu_b, w->wu_t, w->wu, w, k, i);
                    calcFPRes(w->Rwu, w->wu_t, w->wu_b, w->l);
                    scaleArray(w->Rwu, sqrt(w->stgs->rho_x), w->n); /* Scaled FP res in ls */

                    nrmRw_con = calcNorm(w->Rwu, w->l);
                    /* K1 */
                    if (w->stgs->k1 && nrmRw_con <= w->stgs->c1 * nrmR_con_old && w->nrmR_con <= r_safe) { // a bit different than matlab
                        memcpy(w->u, w->wu, w->l * sizeof (scs_float));
                        memcpy(w->u_t, w->wu_t, w->l * sizeof (scs_float));
                        memcpy(w->u_b, w->wu_b, w->l * sizeof (scs_float));
                        memcpy(w->R, w->Rwu, w->l * sizeof (scs_float));
                        // TODO: add computation of sb and kapb
                        w->nrmR_con = nrmRw_con;
                        r_safe += w->stgs->sse; // The power already computed at the beginning of the main loop
                        how = 1;
                        break;
                    }
                    /* K2 */
                    if (K2) {
                        slack = nrmRw_con * nrmRw_con - w->stepsize * innerProd(w->dir, w->Rwu, w->l);
                        rhs = w->stgs->sigma * w->nrmR_con * nrmRw_con;
                        if (slack >= rhs) {
                            stepsize2 = (w->stgs->alpha * (slack / (nrmRw_con * nrmRw_con)));
                            addScaledArray(w->u, w->Rwu, w->l, -stepsize2);
                            how = 2;
                            if (r_safe == INFINITY) {
                                r_safe = w->nrmR_con;
                            }
                            break;
                        }
                    }
                } /* end of line-search */
            }

        }
        if (how == -1) { //means that R didn't change
            scaleArray(w->R, 1.0 / sqrt(w->stgs->rho_x), w->n); /* TODO: FP res should be unscaled here, but then needs to be scaled back?? */
            addScaledArray(w->u, w->R, w->l, -w->stgs->alpha);
        }
        if (how != 1) {
            projectLinSysv2(w->u_t, w->u, w, i);
            projectConesv2(w->u_b, w->u_t, w->u, w, k, i); /* u_bar = proj_C(2u_t - u) */
            // TODO: add calculations for sb and kapb
            calcFPRes(w->R, w->u_t, w->u_b, w->l);
            scaleArray(w->R, sqrt(w->stgs->rho_x), w->n);
            w->nrmR_con = calcNorm(w->R, w->l);
        }
    } /* main for loop */

    /* populate solution vectors (unnormalized) and info */
    getSolution(w, sol, info, &r, i);
    info->solveTime = tocq(&solveTimer);

    if (w->stgs->verbose)
        printFooter(d, k, sol, w, info);
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
    if (!d || !k || !info) {
        scs_printf("ERROR: Missing Data, Cone or Info input\n");
        RETURN SCS_NULL;
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
