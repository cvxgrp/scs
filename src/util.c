#include "util.h"
#include "constants.h"

/* return milli-seconds */
#if (defined NOTIMER)

void tic(timer *t) {
}

scs_float tocq(timer *t) {
    return NAN;
}

#elif(defined _WIN32 || _WIN64 || defined _WINDLL)

void tic(timer *t) {
    QueryPerformanceFrequency(&t->freq);
    QueryPerformanceCounter(&t->tic);
}

scs_float tocq(timer *t) {
    QueryPerformanceCounter(&t->toc);
    return (1e3 * (t->toc.QuadPart - t->tic.QuadPart) /
            (scs_float) t->freq.QuadPart);
}
#elif(defined __APPLE__)

void tic(timer *t) {
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

scs_float tocq(timer *t) {
    uint64_t duration; /* elapsed time in clock cycles*/

    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (scs_float) duration / 1e6;
}
#else

void tic(timer *t) {
    clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

scs_float tocq(timer *t) {
    struct timespec temp;

    clock_gettime(CLOCK_MONOTONIC, &t->toc);

    if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
    } else {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
    }
    return (scs_float) temp.tv_sec * 1e3 + (scs_float) temp.tv_nsec / 1e6;
}
#endif

scs_float toc(timer *t) {
    scs_float time = tocq(t);
    scs_printf("time: %8.4f milli-seconds.\n", time);
    return time;
}

scs_float strtoc(char *str, timer *t) {
    scs_float time = tocq(t);
    scs_printf("%s - time: %8.4f milli-seconds.\n", str, time);
    return time;
}

/* LCOV_EXCL_START */
void printConeData(const Cone *k) {
    scs_int i;
    scs_printf("num zeros = %i\n", (int) k->f);
    scs_printf("num LP = %i\n", (int) k->l);
    scs_printf("num SOCs = %i\n", (int) k->qsize);
    scs_printf("soc array:\n");
    for (i = 0; i < k->qsize; i++) {
        scs_printf("%i\n", (int) k->q[i]);
    }
    scs_printf("num SDCs = %i\n", (int) k->ssize);
    scs_printf("sdc array:\n");
    for (i = 0; i < k->ssize; i++) {
        scs_printf("%i\n", (int) k->s[i]);
    }
    scs_printf("num ep = %i\n", (int) k->ep);
    scs_printf("num ed = %i\n", (int) k->ed);
    scs_printf("num PCs = %i\n", (int) k->psize);
    scs_printf("pow array:\n");
    for (i = 0; i < k->psize; i++) {
        scs_printf("%4f\n", (double) k->p[i]);
    }
}

void printWork(const Work *w) {
    scs_int i, l = w->n + w->m;
    scs_printf("\n u_t is \n");
    for (i = 0; i < l; i++) {
        scs_printf("%f\n", w->u_t[i]);
    }
    scs_printf("\n u is \n");
    for (i = 0; i < l; i++) {
        scs_printf("%f\n", w->u[i]);
    }
    scs_printf("\n v is \n");
    for (i = 0; i < l; i++) {
        scs_printf("%f\n", w->v[i]);
    }
}

void printData(const Data *d) {
    scs_printf("m = %i\n", (int) d->m);
    scs_printf("n = %i\n", (int) d->n);

    scs_printf("max_iters = %i\n", (int) d->stgs->max_iters);
    scs_printf("verbose = %i\n", (int) d->stgs->verbose);
    scs_printf("normalize = %i\n", (int) d->stgs->normalize);
    scs_printf("warmStart = %i\n", (int) d->stgs->warm_start);
    scs_printf("eps = %4f\n", d->stgs->eps);
    scs_printf("alpha = %4f\n", d->stgs->alpha);
    scs_printf("rhoX = %4f\n", d->stgs->rho_x);
    scs_printf("cg_rate = %4f\n", d->stgs->cg_rate);
    scs_printf("scale = %4f\n", d->stgs->scale);
}

void printArray(const scs_float *arr, scs_int n, const char *name) {
    scs_int i, j, k = 0;
    scs_int numOnOneLine = 1;
    scs_printf("\n");
    for (i = 0; i < n / numOnOneLine; ++i) {
        for (j = 0; j < numOnOneLine; ++j) {
            scs_printf("%s[%li] = %4f, ", name, (long) k, arr[k]);
            k++;
        }
        scs_printf("\n");
    }
    for (j = k; j < n; ++j) {
        scs_printf("%s[%li] = %4f, ", name, (long) j, arr[j]);
    }
    scs_printf("\n");
}
/* LCOV_EXCL_STOP */

void freeData(Data *d, Cone *k) {
    if (d!=SCS_NULL) {
        if (d->b!=SCS_NULL)
            scs_free(d->b);
        if (d->c!=SCS_NULL)
            scs_free(d->c);
        if (d->stgs!=SCS_NULL)
            scs_free(d->stgs);
        if (d->A!=SCS_NULL) {
            freeAMatrix(d->A);
        }
        scs_free(d);
    }
    if (k!=SCS_NULL) {
        if (k->q!=SCS_NULL)
            scs_free(k->q);
        if (k->s!=SCS_NULL)
            scs_free(k->s);
        if (k->p!=SCS_NULL)
            scs_free(k->p);
        scs_free(k);
    }
}

void freeSol(Sol *sol) {
    if (sol!=SCS_NULL) {
        if (sol->x!=SCS_NULL) {
            scs_free(sol->x);
        }
        if (sol->y!=SCS_NULL) {
            scs_free(sol->y);
        }
        if (sol->s!=SCS_NULL) {
            scs_free(sol->s);
        }
        scs_free(sol);
    }
}

void freeInfo(Info *info) {
    if (info != SCS_NULL) {
        if (info->progress_iter != SCS_NULL){
            scs_free(info->progress_iter);
        }
        if (info->progress_relgap != SCS_NULL){
            scs_free(info->progress_relgap);
        }
        if (info->progress_resdual != SCS_NULL){
            scs_free(info->progress_resdual);
        }
        if (info->progress_respri != SCS_NULL){
            scs_free(info->progress_respri);
        }
        if (info->progress_pcost != SCS_NULL){
            scs_free(info->progress_pcost);
        }
        if (info->progress_dcost != SCS_NULL){
            scs_free(info->progress_dcost);
        }
        if (info->progress_norm_fpr != SCS_NULL){
            scs_free(info->progress_norm_fpr);
        }
        scs_free(info);
    }
}

/* assumes d->stgs already allocated memory */
void setDefaultSettings(Data *d) {
    d->stgs->max_iters = MAX_ITERS; /* maximum iterations to take: 2500 */
    d->stgs->previous_max_iters = -1; /* maximum iterations of previous invocation */
    d->stgs->eps = EPS; /* convergence tolerance: 1e-3 */
    d->stgs->alpha = ALPHA; /* relaxation parameter: 1.5 */
    d->stgs->rho_x = RHO_X; /* parameter rho_x: 1e-3 */
    d->stgs->scale = SCALE; /* if normalized, rescales by this factor: 1 */
    d->stgs->cg_rate = CG_RATE; /* for indirect, tolerance goes down like (1/iter)^CG_RATE: 2 */
    d->stgs->verbose = VERBOSE; /* boolean, write out progress: 1 */
    d->stgs->normalize = NORMALIZE; /* boolean, heuristic data rescaling: 1 */
    d->stgs->warm_start = WARM_START;
    
    /* -----------------------------
     * SuperSCS-specific parameters
     * ----------------------------- */
    d->stgs->beta = BETA_DEFAULT;
    d->stgs->c1 = C1_DEFAULT;
    d->stgs->c_bl = C_BL_DEFAULT;
    d->stgs->k0 = K0_DEFAULT;
    d->stgs->k1 = K1_DEFAULT;
    d->stgs->k2 = K2_DEFAULT;
    d->stgs->ls = LS_DEFAULT;
    d->stgs->sigma = SIGMA_DEFAULT;
    d->stgs->thetabar = THETABAR_DEFAULT;
    d->stgs->sse = SSE_DEFAULT;
    d->stgs->memory = MEMORY_DEFAULT;
    d->stgs->direction = DIRECTION_DEFAULT;
    d->stgs->do_super_scs = 1; /* whether to run in SuperSCS mode (default: 1) */
    d->stgs->do_record_progress = 0;
}
