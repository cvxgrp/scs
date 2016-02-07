#include "private.h"

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

char *getLinSysMethod(const AMatrix *A, const Settings *s) {
    char *str = scs_malloc(sizeof(char) * 128);
    sprintf(str, "sparse-indirect, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)",
            (long)A->p[A->n], s->cg_rate);
    return str;
}

char *getLinSysSummary(Priv *p, const Info *info) {
    char *str = scs_malloc(sizeof(char) * 128);
    sprintf(str,
            "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
            (scs_float)p->totCgIts / (info->iter + 1),
            p->totalSolveTime / (info->iter + 1) / 1e3);
    p->totCgIts = 0;
    p->totalSolveTime = 0;
    return str;
}

/* M = inv ( diag ( RHO_X * I + A'A ) ) */
void getPreconditioner(const AMatrix *A, const Settings *stgs, Priv *p) {
    scs_int i;
    scs_float *M = p->M;

#if EXTRAVERBOSE > 0
    scs_printf("getting pre-conditioner\n");
#endif

    for (i = 0; i < A->n; ++i) {
        M[i] = 1 / (stgs->rho_x +
                    calcNormSq(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
        /* M[i] = 1; */
    }

#if EXTRAVERBOSE > 0
    scs_printf("finished getting pre-conditioner\n");
#endif
}

static void transpose(const AMatrix *A, Priv *p) {
    scs_int *Ci = p->At->i;
    scs_int *Cp = p->At->p;
    scs_float *Cx = p->At->x;
    scs_int m = A->m;
    scs_int n = A->n;

    scs_int *Ap = A->p;
    scs_int *Ai = A->i;
    scs_float *Ax = A->x;

    scs_int i, j, q, *z, c1, c2;
#if EXTRAVERBOSE > 0
    timer transposeTimer;
    scs_printf("transposing A\n");
    tic(&transposeTimer);
#endif

    z = scs_calloc(m, sizeof(scs_int));
    for (i = 0; i < Ap[n]; i++)
        z[Ai[i]]++;      /* row counts */
    cs_cumsum(Cp, z, m); /* row pointers */

    for (j = 0; j < n; j++) {
        c1 = Ap[j];
        c2 = Ap[j + 1];
        for (i = c1; i < c2; i++) {
            q = z[Ai[i]];
            Ci[q] = j; /* place A(i,j) as entry C(j,i) */
            Cx[q] = Ax[i];
            z[Ai[i]]++;
        }
    }
    scs_free(z);

#if EXTRAVERBOSE > 0
    scs_printf("finished transposing A, time: %1.2es\n",
               tocq(&transposeTimer) / 1e3);
#endif
}

void freePriv(Priv *p) {
    if (p) {
        if (p->p)
            scs_free(p->p);
        if (p->r)
            scs_free(p->r);
        if (p->Gp)
            scs_free(p->Gp);
        if (p->tmp)
            scs_free(p->tmp);
        if (p->At) {
            AMatrix *At = p->At;
            if (At->i)
                scs_free(At->i);
            if (At->x)
                scs_free(At->x);
            if (At->p)
                scs_free(At->p);
            scs_free(At);
        }
        if (p->z)
            scs_free(p->z);
        if (p->M)
            scs_free(p->M);
        scs_free(p);
    }
}

/*y = (RHO_X * I + A'A)x */
static void matVec(const AMatrix *A, const Settings *s, Priv *p,
                   const scs_float *x, scs_float *y) {
    scs_float *tmp = p->tmp;
    memset(tmp, 0, A->m * sizeof(scs_float));
    accumByA(A, p, x, tmp);
    memset(y, 0, A->n * sizeof(scs_float));
    accumByAtrans(A, p, tmp, y);
    addScaledArray(y, x, A->n, s->rho_x);
}

void accumByAtrans(const AMatrix *A, Priv *p, const scs_float *x,
                   scs_float *y) {
    _accumByAtrans(A->n, A->x, A->i, A->p, x, y);
}

void accumByA(const AMatrix *A, Priv *p, const scs_float *x, scs_float *y) {
    _accumByAtrans(p->At->n, p->At->x, p->At->i, p->At->p, x, y);
}

static void applyPreConditioner(scs_float *M, scs_float *z, scs_float *r,
                                scs_int n, scs_float *ipzr) {
    scs_int i;
    *ipzr = 0;
    for (i = 0; i < n; ++i) {
        z[i] = r[i] * M[i];
        *ipzr += z[i] * r[i];
    }
}

Priv *initPriv(const AMatrix *A, const Settings *stgs) {
    Priv *p = scs_calloc(1, sizeof(Priv));
    p->p = scs_malloc((A->n) * sizeof(scs_float));
    p->r = scs_malloc((A->n) * sizeof(scs_float));
    p->Gp = scs_malloc((A->n) * sizeof(scs_float));
    p->tmp = scs_malloc((A->m) * sizeof(scs_float));

    /* memory for A transpose */
    p->At = scs_malloc(sizeof(AMatrix));
    p->At->m = A->n;
    p->At->n = A->m;
    p->At->i = scs_malloc((A->p[A->n]) * sizeof(scs_int));
    p->At->p = scs_malloc((A->m + 1) * sizeof(scs_int));
    p->At->x = scs_malloc((A->p[A->n]) * sizeof(scs_float));
    transpose(A, p);

    /* preconditioner memory */
    p->z = scs_malloc((A->n) * sizeof(scs_float));
    p->M = scs_malloc((A->n) * sizeof(scs_float));
    getPreconditioner(A, stgs, p);

    p->totalSolveTime = 0;
    p->totCgIts = 0;
    if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At || !p->At->i ||
        !p->At->p || !p->At->x) {
        freePriv(p);
        return SCS_NULL;
    }
    return p;
}

/* solves (I+A'A)x = b, s warm start, solution stored in b */
static scs_int pcg(const AMatrix *A, const Settings *stgs, Priv *pr,
                   const scs_float *s, scs_float *b, scs_int max_its,
                   scs_float tol) {
    scs_int i, n = A->n;
    scs_float ipzr, ipzrOld, alpha;
    scs_float *p = pr->p;   /* cg direction */
    scs_float *Gp = pr->Gp; /* updated CG direction */
    scs_float *r = pr->r;   /* cg residual */
    scs_float *z = pr->z;   /* for preconditioning */
    scs_float *M = pr->M;   /* inverse diagonal preconditioner */

    if (s == SCS_NULL) {
        memcpy(r, b, n * sizeof(scs_float));
        memset(b, 0, n * sizeof(scs_float));
    } else {
        matVec(A, stgs, pr, s, r);
        addScaledArray(r, b, n, -1);
        scaleArray(r, -1, n);
        memcpy(b, s, n * sizeof(scs_float));
    }

    /* check to see if we need to run CG at all */
    if (calcNorm(r, n) < MIN(tol, 1e-18)) {
        return 0;
    }

    applyPreConditioner(M, z, r, n, &ipzr);
    memcpy(p, z, n * sizeof(scs_float));

    for (i = 0; i < max_its; ++i) {
        matVec(A, stgs, pr, p, Gp);
        alpha = ipzr / innerProd(p, Gp, n);
        addScaledArray(b, p, n, alpha);
        addScaledArray(r, Gp, n, -alpha);

        if (calcNorm(r, n) < tol) {
#if EXTRAVERBOSE > 0
            scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol,
                       calcNorm(r, n), (long)i + 1);
#endif
            return i + 1;
        }
        ipzrOld = ipzr;
        applyPreConditioner(M, z, r, n, &ipzr);

        scaleArray(p, ipzr / ipzrOld, n);
        addScaledArray(p, z, n, 1);
    }
    return i;
}

scs_int solveLinSys(const AMatrix *A, const Settings *stgs, Priv *p,
                    scs_float *b, const scs_float *s, scs_int iter) {
    scs_int cgIts;
    timer linsysTimer;
    scs_float cgTol =
        calcNorm(b, A->n) *
        (iter < 0 ? CG_BEST_TOL
                  : CG_MIN_TOL / POWF((scs_float)iter + 1, stgs->cg_rate));

    tic(&linsysTimer);
    /* solves Mx = b, for x but stores result in b */
    /* s contains warm-start (if available) */
    accumByAtrans(A, p, &(b[A->n]), b);
    /* solves (I+A'A)x = b, s warm start, solution stored in b */
    cgIts = pcg(A, stgs, p, s, b, A->n, MAX(cgTol, CG_BEST_TOL));
    scaleArray(&(b[A->n]), -1, A->m);
    accumByA(A, p, b, &(b[A->n]));

    if (iter >= 0) {
        p->totCgIts += cgIts;
    }

    p->totalSolveTime += tocq(&linsysTimer);
#if EXTRAVERBOSE > 0
    scs_printf("linsys solve time: %1.2es\n", tocq(&linsysTimer) / 1e3);
#endif
    return 0;
}
