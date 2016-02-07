#include "scs.h"
#include "normalize.h"

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
    " Iter ",    " pri res ", " dua res ", " rel gap ",
    " pri obj ", " dua obj ", " kap/tau ", " time (s)",
};
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 8;
static const scs_int LINE_LEN = 76;

static scs_int scs_isnan(scs_float x) {
    DEBUG_FUNC
    RETURN(x == NAN || x != x);
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
                   stgs->eps, stgs->alpha, (int)stgs->max_iters,
                   (int)stgs->normalize, stgs->scale);
    } else {
        scs_printf("eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i\n",
                   stgs->eps, stgs->alpha, (int)stgs->max_iters,
                   (int)stgs->normalize);
    }
    scs_printf("Variables n = %i, constraints m = %i\n", (int)d->n, (int)d->m);
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
                sol->x = scs_malloc(sizeof(scs_float) * n);
            scaleArray(sol->x, NAN, n);
        }
        if (m > 0) {
            if (!sol->y)
                sol->y = scs_malloc(sizeof(scs_float) * m);
            scaleArray(sol->y, NAN, m);
            if (!sol->s)
                sol->s = scs_malloc(sizeof(scs_float) * m);
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
    memset(w->v, 0, n * sizeof(scs_float));
    memcpy(w->u, sol->x, n * sizeof(scs_float));
    memcpy(&(w->u[n]), sol->y, m * sizeof(scs_float));
    memcpy(&(w->v[n]), sol->s, m * sizeof(scs_float));
    w->u[n + m] = 1.0;
    w->v[n + m] = 0.0;
#ifndef NOVALIDATE
    for (i = 0; i < n + m + 1; ++i) {
        if (scs_isnan(w->u[i]))
            w->u[i] = 0;
        if (scs_isnan(w->v[i]))
            w->v[i] = 0;
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
    memset(pr, 0, w->m * sizeof(scs_float));
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
    memset(dr, 0, w->n * sizeof(scs_float));
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
    scs_int l = w->n + w->m + 1;
    memset(w->u, 0, l * sizeof(scs_float));
    memset(w->v, 0, l * sizeof(scs_float));
    w->u[l - 1] = SQRTF((scs_float)l);
    w->v[l - 1] = SQRTF((scs_float)l);
    RETURN;
}

/* status < 0 indicates failure */
static scs_int projectLinSys(Work *w, scs_int iter) {
    /* ut = u + v */
    DEBUG_FUNC
    scs_int n = w->n, m = w->m, l = n + m + 1, status;
    memcpy(w->u_t, w->u, l * sizeof(scs_float));
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

void printSol(Work *w, Sol *sol, Info *info) {
    DEBUG_FUNC
    scs_int i;
    scs_printf("%s\n", info->status);
    if (sol->x != SCS_NULL) {
        for (i = 0; i < w->n; ++i) {
            scs_printf("x[%i] = %4f\n", (int)i, sol->x[i]);
        }
    }
    if (sol->y != SCS_NULL) {
        for (i = 0; i < w->m; ++i) {
            scs_printf("y[%i] = %4f\n", (int)i, sol->y[i]);
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
        sol->y = scs_malloc(sizeof(scs_float) * w->m);
    }
    memcpy(sol->y, &(w->u[w->n]), w->m * sizeof(scs_float));
    RETURN;
}

static void sets(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->s) {
        sol->s = scs_malloc(sizeof(scs_float) * w->m);
    }
    memcpy(sol->s, &(w->v[w->n]), w->m * sizeof(scs_float));
    RETURN;
}

static void setx(Work *w, Sol *sol) {
    DEBUG_FUNC
    if (!sol->x)
        sol->x = scs_malloc(sizeof(scs_float) * w->n);
    memcpy(sol->x, w->u, w->n * sizeof(scs_float));
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
                   INDETERMINATE_TOL * SQRTF((scs_float)l)) {
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
    scs_printf("%*i|", (int)strlen(HEADER[0]), (int)i);
    scs_printf("%*.2e ", (int)HSPACE, r->resPri);
    scs_printf("%*.2e ", (int)HSPACE, r->resDual);
    scs_printf("%*.2e ", (int)HSPACE, r->relGap);
    scs_printf("%*.2e ", (int)HSPACE, r->cTx_by_tau / r->tau);
    scs_printf("%*.2e ", (int)HSPACE, -r->bTy_by_tau / r->tau);
    scs_printf("%*.2e ", (int)HSPACE, r->kap / r->tau);
    scs_printf("%*.2e ", (int)HSPACE, tocq(solveTimer) / 1e3);
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
    scs_float *t = scs_malloc(sizeof(scs_float) * m);
    memcpy(t, y, m * sizeof(scs_float));
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
    scs_float *t = scs_malloc(sizeof(scs_float) * m);
    memcpy(t, s, m * sizeof(scs_float));
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
                   (long)d->m, (long)d->n);
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
        scs_printf("max_iters must be positive\n");
        RETURN - 1;
    }
    if (stgs->eps <= 0) {
        scs_printf("eps tolerance must be positive\n");
        RETURN - 1;
    }
    if (stgs->alpha <= 0 || stgs->alpha >= 2) {
        scs_printf("alpha must be in (0,2)\n");
        RETURN - 1;
    }
    if (stgs->rho_x <= 0) {
        scs_printf("rhoX must be positive (1e-3 works well).\n");
        RETURN - 1;
    }
    if (stgs->scale <= 0) {
        scs_printf("scale must be positive (1 works well).\n");
        RETURN - 1;
    }
    RETURN 0;
}

static Work *initWork(const Data *d, const Cone *k) {
    DEBUG_FUNC
    Work *w = scs_calloc(1, sizeof(Work));
    scs_int l = d->n + d->m + 1;
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
    /* allocate workspace: */
    w->u = scs_malloc(l * sizeof(scs_float));
    w->v = scs_malloc(l * sizeof(scs_float));
    w->u_t = scs_malloc(l * sizeof(scs_float));
    w->u_prev = scs_malloc(l * sizeof(scs_float));
    w->h = scs_malloc((l - 1) * sizeof(scs_float));
    w->g = scs_malloc((l - 1) * sizeof(scs_float));
    w->pr = scs_malloc(d->m * sizeof(scs_float));
    w->dr = scs_malloc(d->n * sizeof(scs_float));
    w->b = scs_malloc(d->m * sizeof(scs_float));
    w->c = scs_malloc(d->n * sizeof(scs_float));
    if (!w->u || !w->v || !w->u_t || !w->u_prev || !w->h || !w->g || !w->pr ||
        !w->dr || !w->b || !w->c) {
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
        w->scal = scs_malloc(sizeof(Scaling));
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
    memcpy(w->b, d->b, d->m * sizeof(scs_float));
    memcpy(w->c, d->c, d->n * sizeof(scs_float));

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
    memcpy(w->h, w->c, n * sizeof(scs_float));
    memcpy(&(w->h[n]), w->b, m * sizeof(scs_float));
    memcpy(w->g, w->h, (n + m) * sizeof(scs_float));
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
        memcpy(w->u_prev, w->u, (w->n + w->m + 1) * sizeof(scs_float));

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

void scs_finish(Work *w) {
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

Work *scs_init(const Data *d, const Cone *k, Info *info) {
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
scs_int scs(const Data *d, const Cone *k, Sol *sol, Info *info) {
    DEBUG_FUNC
    scs_int status;
#if (defined _WIN32 || defined _WIN64) && !defined MATLAB_MEX_FILE &&          \
    !defined PYTHON
    /* sets width of exponent for floating point numbers to 2 instead of 3 */
    unsigned int old_output_format = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    Work *w = scs_init(d, k, info);
#if EXTRAVERBOSE > 0
    scs_printf("size of scs_int = %lu, size of scs_float = %lu\n",
               sizeof(scs_int), sizeof(scs_float));
#endif
    if (w) {
        scs_solve(w, d, k, sol, info);
        status = info->statusVal;
    } else {
        status = failure(SCS_NULL, d ? d->m : -1, d ? d->n : -1, sol, info,
                         SCS_FAILED, "could not initialize work", "Failure");
    }
    scs_finish(w);
    RETURN status;
}
