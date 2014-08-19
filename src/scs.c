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

/* printing header */
static const char* HEADER[] = { " Iter ", " pri res ", " dua res ", " rel gap ", " pri obj ", " dua obj ", " kap/tau ",
		" time (s)", };
static const scs_int HSPACE = 9;
static const scs_int HEADER_LEN = 8;
static scs_int _lineLen_;

static scs_int scs_isnan(scs_float x) {
	return (x == NAN || x != x);
}

static void printInitHeader(Data * d, Work * w, Cone * k) {
	scs_int i;
	char * coneStr = getConeHeader(k);
	char * linSysMethod = getLinSysMethod(d, w->p);
	_lineLen_ = -1;
	for (i = 0; i < HEADER_LEN; ++i) {
		_lineLen_ += (scs_int) strlen(HEADER[i]) + 1;
	}
	for (i = 0; i < _lineLen_; ++i) {
		scs_printf("-");
	}
	scs_printf("\n\tSCS v%s - Splitting Conic Solver\n\t(c) Brendan O'Donoghue, Stanford University, 2012\n",
			SCS_VERSION);
	for (i = 0; i < _lineLen_; ++i) {
		scs_printf("-");
	}
    scs_printf("\n");
	if (linSysMethod) {
		scs_printf("Lin-sys: %s\n", linSysMethod);
		scs_free(linSysMethod);
	}
	if (d->normalize) {
		scs_printf("eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i, scale = %2.2f\n", d->eps, d->alpha,
				(int) d->max_iters, (int) d->normalize, d->scale);
	} else {
		scs_printf("eps = %.2e, alpha = %.2f, max_iters = %i, normalize = %i\n", d->eps, d->alpha, (int) d->max_iters,
				(int) d->normalize);
	}
	scs_printf("Variables n = %i, constraints m = %i\n", (int) d->n, (int) d->m);
	scs_printf("%s", coneStr);
	scs_free(coneStr);
#ifdef MATLAB_MEX_FILE
	mexEvalString("drawnow;");
#endif
}

static scs_int failureDefaultReturn(Data * d, Work * w, Sol * sol, Info * info, char * msg) {
	info->relGap = NAN;
	info->resPri = NAN;
	info->resDual = NAN;
	info->pobj = NAN;
	info->dobj = NAN;
	info->iter = -1;
	info->statusVal = FAILURE;
	info->solveTime = NAN;
	strcpy(info->status, "Failure");
	if (!sol->x)
		sol->x = scs_malloc(sizeof(scs_float) * d->n);
	scaleArray(sol->x, NAN, d->n);
	if (!sol->y)
		sol->y = scs_malloc(sizeof(scs_float) * d->m);
	scaleArray(sol->y, NAN, d->m);
	if (!sol->s)
		sol->s = scs_malloc(sizeof(scs_float) * d->m);
	scaleArray(sol->s, NAN, d->m);
	scs_printf("FAILURE:%s\n", msg);
	return FAILURE;
}

static void warmStartVars(Data * d, Work * w, Sol * sol) {
	scs_int i, n = d->n, m = d->m;
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
	if (d->normalize)
		normalizeWarmStart(d, w);
}

static scs_float calcPrimalResid(Data * d, Work * w, scs_float * x, scs_float * s, scs_float tau, scs_float *nmAxs) {
	scs_int i;
	scs_float pres = 0, scale, *pr = w->pr, *D = w->D;
	*nmAxs = 0;
	memset(pr, 0, d->m * sizeof(scs_float));
	accumByA(d, w->p, x, pr);
	addScaledArray(pr, s, d->m, 1.0); /* pr = Ax + s */
	for (i = 0; i < d->m; ++i) {
		scale = d->normalize ? D[i] / (w->sc_b * d->scale) : 1;
		scale = scale * scale;
		*nmAxs += (pr[i] * pr[i]) * scale;
		pres += (pr[i] - d->b[i] * tau) * (pr[i] - d->b[i] * tau) * scale;
	}
	*nmAxs = SQRTF(*nmAxs);
	return SQRTF(pres); /* norm(Ax + s - b * tau) */
}

static scs_float calcDualResid(Data * d, Work * w, scs_float * y, scs_float tau, scs_float *nmATy) {
	scs_int i;
	scs_float dres = 0, scale, *dr = w->dr, *E = w->E;
	*nmATy = 0;
	memset(dr, 0, d->n * sizeof(scs_float));
	accumByAtrans(d, w->p, y, dr); /* dr = A'y */
	for (i = 0; i < d->n; ++i) {
		scale = d->normalize ? E[i] / (w->sc_c * d->scale) : 1;
		scale = scale * scale;
		*nmATy += (dr[i] * dr[i]) * scale;
		dres += (dr[i] + d->c[i] * tau) * (dr[i] + d->c[i] * tau) * scale;
	}
	*nmATy = SQRTF(*nmATy);
	return SQRTF(dres); /* norm(A'y + c * tau) */
}

static scs_float fastCalcPrimalResid(Data * d, Work * w, scs_float * nmAxs) {
	scs_int i, n = d->n, m = d->m;
	scs_float pres = 0, scale, *pr = w->pr, *D = w->D, tau = ABS(w->u[n + m]);
	*nmAxs = 0;
	memcpy(pr, &(w->u[n]), m * sizeof(scs_float)); /* overwrite pr */
	addScaledArray(pr, &(w->u_prev[n]), m, d->alpha - 2);
	addScaledArray(pr, &(w->u_t[n]), m, 1 - d->alpha);
	addScaledArray(pr, d->b, m, w->u_t[n + m]); /* pr = Ax + s */
	for (i = 0; i < m; ++i) {
		scale = d->normalize ? D[i] / (w->sc_b * d->scale) : 1;
		scale = scale * scale;
		*nmAxs += (pr[i] * pr[i]) * scale;
		pres += (pr[i] - d->b[i] * tau) * (pr[i] - d->b[i] * tau) * scale;
	}
	*nmAxs = SQRTF(*nmAxs);
	return SQRTF(pres); /* norm(Ax + s - b * tau) */
}

static void getInfo(Data * d, Work * w, Sol * sol, Info * info) {
	scs_float cTx, bTy, nmAxs, nmATy, nmpr, nmdr;
	scs_float * x = sol->x, *y = sol->y, *s = sol->s;

	/* unNomalized */
	nmpr = calcPrimalResid(d, w, x, s, 1, &nmAxs); /* pr = Ax + s - b */
	nmdr = calcDualResid(d, w, y, 1, &nmATy); /* dr = A'y + c */

	cTx = innerProd(x, d->c, d->n);
	bTy = innerProd(y, d->b, d->m);
	if (d->normalize) {
		cTx /= (d->scale * w->sc_c * w->sc_b);
		bTy /= (d->scale * w->sc_c * w->sc_b);
	}
	info->pobj = cTx;
	info->dobj = -bTy;
	if (info->statusVal == SOLVED) {
		info->relGap = ABS(cTx + bTy) / (1 + ABS(cTx) + ABS(bTy));
		info->resPri = nmpr / (1 + w->nm_b);
		info->resDual = nmdr / (1 + w->nm_c);
	} else {
		if (info->statusVal == UNBOUNDED) {
			info->dobj = NAN;
			info->relGap = NAN;
			info->resPri = w->nm_c * nmAxs / -cTx;
			info->resDual = NAN;
			scaleArray(x, -1 / cTx, d->n);
			scaleArray(s, -1 / cTx, d->m);
			info->pobj = -1;
		} else {
			info->pobj = NAN;
			info->relGap = NAN;
			info->resPri = NAN;
			info->resDual = w->nm_b * nmATy / -bTy;
			scaleArray(y, -1 / bTy, d->m);
			info->dobj = -1;
		}
	}
}

static void coldStartVars(Data * d, Work * w) {
	scs_int l = d->n + d->m + 1;
	memset(w->u, 0, l * sizeof(scs_float));
	memset(w->v, 0, l * sizeof(scs_float));
	w->u[l - 1] = SQRTF((scs_float) l);
	w->v[l - 1] = SQRTF((scs_float) l);
}

/* status < 0 indicates failure */
static scs_int projectLinSys(Data * d, Work * w, scs_int iter) {
	/* ut = u + v */
	scs_int n = d->n, m = d->m, l = n + m + 1, status;
	memcpy(w->u_t, w->u, l * sizeof(scs_float));
	addScaledArray(w->u_t, w->v, l, 1.0);

	scaleArray(w->u_t, d->rho_x, n);

	addScaledArray(w->u_t, w->h, l - 1, -w->u_t[l - 1]);
	addScaledArray(w->u_t, w->h, l - 1, -innerProd(w->u_t, w->g, l - 1) / (w->gTh + 1));
	scaleArray(&(w->u_t[n]), -1, m);

	status = solveLinSys(d, w->p, w->u_t, w->u, iter);

	w->u_t[l - 1] += innerProd(w->u_t, w->h, l - 1);

	return status;
}

static void freeWork(Work * w) {
	if (w) {
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
		if (w->D)
			scs_free(w->D);
		if (w->E)
			scs_free(w->E);
		if (w->pr)
			scs_free(w->pr);
		if (w->dr)
			scs_free(w->dr);
		scs_free(w);
	}
}

void printSol(Data * d, Sol * sol, Info * info) {
	scs_int i;
	scs_printf("%s\n", info->status);
	if (sol->x != NULL) {
		for (i = 0; i < d->n; ++i) {
			scs_printf("x[%i] = %4f\n", (int) i, sol->x[i]);
		}
	}
	if (sol->y != NULL) {
		for (i = 0; i < d->m; ++i) {
			scs_printf("y[%i] = %4f\n", (int) i, sol->y[i]);
		}
	}
}

static void updateDualVars(Data * d, Work * w) {
	scs_int i, n = d->n, l = n + d->m + 1;
	/* this does not relax 'x' variable */
	for (i = n; i < l; ++i) {
		w->v[i] += (w->u[i] - d->alpha * w->u_t[i] - (1.0 - d->alpha) * w->u_prev[i]);
	}
}

/* status < 0 indicates failure */
static scs_int projectCones(Data *d, Work * w, Cone * k, scs_int iter) {
	scs_int i, n = d->n, l = n + d->m + 1, status;
	/* this does not relax 'x' variable */
	for (i = 0; i < n; ++i) {
		w->u[i] = w->u_t[i] - w->v[i];
	}
	for (i = n; i < l; ++i) {
		w->u[i] = d->alpha * w->u_t[i] + (1 - d->alpha) * w->u_prev[i] - w->v[i];
	}
	/* u = [x;y;tau] */
	status = projDualCone(&(w->u[n]), k, &(w->u_prev[n]), iter);
	if (w->u[l - 1] < 0.0)
		w->u[l - 1] = 0.0;

	return status;
}

static scs_int solved(Data * d, Sol * sol, Info * info, scs_float tau) {
	strcpy(info->status, info->statusVal != 0 ? "Solved" : "Solved/Inaccurate");
	scaleArray(sol->x, 1.0 / tau, d->n);
	scaleArray(sol->y, 1.0 / tau, d->m);
	scaleArray(sol->s, 1.0 / tau, d->m);
	return SOLVED;
}

static scs_int indeterminate(Data * d, Sol * sol, Info * info) {
	strcpy(info->status, "Indeterminate");
	scaleArray(sol->x, NAN, d->n);
	scaleArray(sol->y, NAN, d->m);
	scaleArray(sol->s, NAN, d->m);
	return INDETERMINATE;
}

static scs_int infeasible(Data * d, Sol * sol, Info * info) {
	strcpy(info->status, info->statusVal != 0 ? "Infeasible" : "Infeasible/Inaccurate");
	/*scaleArray(sol->y,-1/ip_y,d->m); */
	scaleArray(sol->x, NAN, d->n);
	scaleArray(sol->s, NAN, d->m);
	return INFEASIBLE;
}

static scs_int unbounded(Data * d, Sol * sol, Info * info) {
	strcpy(info->status, info->statusVal != 0 ? "Unbounded" : "Unbounded/Inaccurate");
	/*scaleArray(sol->x,-1/ip_x,d->n); */
	scaleArray(sol->y, NAN, d->m);
	return UNBOUNDED;
}

static void sety(Data * d, Work * w, Sol * sol) {
	if (!sol->y)
		sol->y = scs_malloc(sizeof(scs_float) * d->m);
	memcpy(sol->y, &(w->u[d->n]), d->m * sizeof(scs_float));
}

static void sets(Data * d, Work * w, Sol * sol) {
	if (!sol->s)
		sol->s = scs_malloc(sizeof(scs_float) * d->m);
	memcpy(sol->s, &(w->v[d->n]), d->m * sizeof(scs_float));
}

static void setx(Data * d, Work * w, Sol * sol) {
	if (!sol->x)
		sol->x = scs_malloc(sizeof(scs_float) * d->n);
	memcpy(sol->x, w->u, d->n * sizeof(scs_float));
}

static void setSolution(Data * d, Work * w, Sol * sol, Info * info) {
	scs_int l = d->n + d->m + 1;
	setx(d, w, sol);
	sety(d, w, sol);
	sets(d, w, sol);
	if (info->statusVal == 0 || info->statusVal == SOLVED) {
		scs_float tau = w->u[l - 1];
		scs_float kap = ABS(w->v[l - 1]);
		if (tau > INDETERMINATE_TOL && tau > kap) {
			info->statusVal = solved(d, sol, info, tau);
		} else {
			if (calcNorm(w->u, l) < INDETERMINATE_TOL * SQRTF((scs_float) l)) {
				info->statusVal = indeterminate(d, sol, info);
			} else {
				scs_float bTy = innerProd(d->b, sol->y, d->m);
				scs_float cTx = innerProd(d->c, sol->x, d->n);
				if (bTy < cTx) {
					info->statusVal = infeasible(d, sol, info);
				} else {
					info->statusVal = unbounded(d, sol, info);
				}
			}
		}
	} else if (info->statusVal == INFEASIBLE) {
		info->statusVal = infeasible(d, sol, info);
	} else {
		info->statusVal = unbounded(d, sol, info);
	}
}

static void printSummary(scs_int i, struct residuals *r, timer * solveTimer) {
	scs_printf("%*i|", (int) strlen(HEADER[0]), (int) i);
	scs_printf("%*.2e ", (int) HSPACE, r->resPri);
	scs_printf("%*.2e ", (int) HSPACE, r->resDual);
	scs_printf("%*.2e ", (int) HSPACE, r->relGap);
	scs_printf("%*.2e ", (int) HSPACE, r->cTx);
	scs_printf("%*.2e ", (int) HSPACE, -r->bTy);
	scs_printf("%*.2e ", (int) HSPACE, r->kap / r->tau);
	scs_printf("%*.2e ", (int) HSPACE, tocq(solveTimer) / 1e3);
	scs_printf("\n");
#ifdef MATLAB_MEX_FILE
	mexEvalString("drawnow;");
#endif
}

static void printHeader(Data * d, Work * w, Cone * k) {
	scs_int i;
    if (d->warm_start)
        scs_printf("SCS using variable warm-starting\n");
    for (i = 0; i < _lineLen_; ++i) {
        scs_printf("-");
    }
    scs_printf("\n");
    for (i = 0; i < HEADER_LEN - 1; ++i) {
        scs_printf("%s|", HEADER[i]);
    }
    scs_printf("%s\n", HEADER[HEADER_LEN - 1]);
    for (i = 0; i < _lineLen_; ++i) {
        scs_printf("-");
    }
    scs_printf("\n");
#ifdef MATLAB_MEX_FILE
    mexEvalString("drawnow;");
#endif
}

static void printFooter(Data * d, Work * w, Info * info) {
	scs_int i;
	char * linSysStr = getLinSysSummary(w->p, info);
	char * coneStr = getConeSummary(info);
	for (i = 0; i < _lineLen_; ++i) {
		scs_printf("-");
	}
	scs_printf("\nStatus: %s\n", info->status);
	if (info->iter == d->max_iters) {
		scs_printf("Hit max_iters, solution may be inaccurate\n");
	}
	scs_printf("Timing: Total solve time: %1.2es\n", info->solveTime / 1e3);

	if (linSysStr) {
		scs_printf("%s", linSysStr);
		scs_free(linSysStr);
	}
	if (coneStr) {
		scs_printf("%s", coneStr);
		scs_free(coneStr);
	}

	for (i = 0; i < _lineLen_; ++i) {
		scs_printf("-");
	}
	scs_printf("\n");

	if (info->statusVal == INFEASIBLE) {
		scs_printf("Certificate of primal infeasibility:\n");
		scs_printf("|A'y|_2 * |b|_2 = %.4e\n", info->resDual);
		scs_printf("dist(y, K*) = 0\n");
		scs_printf("b'y = %.4f\n", info->dobj);
	} else if (info->statusVal == UNBOUNDED) {
		scs_printf("Certificate of dual infeasibility:\n");
		scs_printf("|Ax + s|_2 * |c|_2 = %.4e\n", info->resPri);
		scs_printf("dist(s, K) = 0\n");
		scs_printf("c'x = %.4f\n", info->pobj);
	} else {
		scs_printf("Error metrics:\n");
		scs_printf("|Ax + s - b|_2 / (1 + |b|_2) = %.4e\n", info->resPri);
		scs_printf("|A'y + c|_2 / (1 + |c|_2) = %.4e\n", info->resDual);
		scs_printf("|c'x + b'y| / (1 + |c'x| + |b'y|) = %.4e\n", info->relGap);
		scs_printf("dist(s, K) = 0, dist(y, K*) = 0, s'y = 0\n");
		for (i = 0; i < _lineLen_; ++i) {
			scs_printf("-");
		}
		scs_printf("\n");
		scs_printf("c'x = %.4f, -b'y = %.4f\n", info->pobj, info->dobj);
	}
	for (i = 0; i < _lineLen_; ++i) {
		scs_printf("=");
	}
	scs_printf("\n");
#ifdef MATLAB_MEX_FILE
	mexEvalString("drawnow;");
#endif
}

static scs_int converged(Data * d, Work * w, struct residuals * r, scs_int iter) {
	scs_float nmpr, nmdr, tau, kap, *x, *y, cTx, nmAxs, bTy, nmATy, rpri, rdua, gap;
	scs_int n = d->n, m = d->m;
	if (iter % CONVERGED_INTERVAL != 0) {
		return 0;
	}
	x = w->u;
	y = &(w->u[n]);
	tau = ABS(w->u[n + m]);
	kap = ABS(w->v[n + m]) / (d->normalize ? (d->scale * w->sc_c * w->sc_b) : 1);
    r->tau = tau;
	r->kap = kap;

	/* requires mult by A:
	 nmpr = calcPrimalResid(d, w, w->u, &(w->v[n]), ABS(w->u[n + m]), &nmAxs);
	 */

	/* does not require mult by A: */
	nmpr = fastCalcPrimalResid(d, w, &nmAxs);
	cTx = innerProd(x, d->c, n) / (d->normalize ? (d->scale * w->sc_c * w->sc_b) : 1);

	r->resPri = cTx < 0 ? w->nm_c * nmAxs / -cTx : NAN;
	if (r->resPri < d->eps) {
		return UNBOUNDED;
	}

	nmdr = calcDualResid(d, w, y, tau, &nmATy);
	bTy = innerProd(y, d->b, m) / (d->normalize ? (d->scale * w->sc_c * w->sc_b) : 1);

	r->resDual = bTy < 0 ? w->nm_b * nmATy / -bTy : NAN;
	if (r->resDual < d->eps) {
		return INFEASIBLE;
	}

	r->cTx = cTx / tau;
	r->bTy = bTy / tau;
	r->relGap = NAN;

	rpri = nmpr / (1 + w->nm_b) / tau;
	rdua = nmdr / (1 + w->nm_c) / tau;
	gap = ABS(cTx + bTy) / (tau + ABS(cTx) + ABS(bTy));
	if (tau > kap) {
		r->resPri = rpri;
		r->resDual = rdua;
		r->relGap = gap;
	}
	return (MAX(MAX(rpri,rdua),gap) < d->eps ? SOLVED : 0);
}

static scs_int validate(Data * d, Cone * k) {
	if (d->m <= 0 || d->n <= 0) {
		scs_printf("m and n must both be greater than 0\n");
		return -1;
	}
    /*
	if (d->m < d->n) {
		scs_printf("m must be greater than or equal to n\n");
		return -1;
	}
    */
	if (validateLinSys(d) < 0) {
		scs_printf("invalid linear system input data\n");
		return -1;
	}
	if (validateCones(d, k) < 0) {
		scs_printf("invalid cone dimensions\n");
		return -1;
	}
	if (d->max_iters <= 0) {
		scs_printf("max_iters must be positive\n");
		return -1;
	}
	if (d->eps <= 0) {
		scs_printf("eps tolerance must be positive\n");
		return -1;
	}
	if (d->alpha <= 0 || d->alpha >= 2) {
		scs_printf("alpha must be in (0,2)\n");
		return -1;
	}
	if (d->rho_x <= 0) {
		scs_printf("rhoX must be positive (1e-3 works well).\n");
		return -1;
	}
	if (d->scale <= 0) {
		scs_printf("scale must be positive (1 works well).\n");
		return -1;
	}
	return 0;
}

static Work * initWork(Data *d, Cone * k) {
	Work * w = scs_calloc(1, sizeof(Work));
	scs_int l = d->n + d->m + 1;
	if (d->verbose) {
		printInitHeader(d, w, k);
	}
	if (!w) {
		scs_printf("ERROR: allocating work failure\n");
		return NULL;
	}
	/* allocate workspace: */
	w->u = scs_malloc(l * sizeof(scs_float));
	w->v = scs_malloc(l * sizeof(scs_float));
	w->u_t = scs_malloc(l * sizeof(scs_float));
	w->u_prev = scs_malloc(l * sizeof(scs_float));
	w->h = scs_malloc((l - 1) * sizeof(scs_float));
	w->g = scs_malloc((l - 1) * sizeof(scs_float));
	w->pr = scs_malloc(d->m * sizeof(scs_float));
	w->dr = scs_malloc(d->n * sizeof(scs_float));
	if (!w->u || !w->v || !w->u_t || !w->u_prev || !w->h || !w->g || !w->pr || !w->dr) {
		scs_printf("ERROR: work memory allocation failure\n");
		scs_finish(d, w);
		return NULL;
	}
	if (d->normalize) {
		normalizeA(d, w, k);
#ifdef EXTRAVERBOSE
	printArray(w->D, d->m, "D");
	scs_printf("norm D = %4f\n", calcNorm(w->D, d->m));
	printArray(w->E, d->n, "E");
	scs_printf("norm E = %4f\n", calcNorm(w->E, d->n));
#endif
	} else {
		w->D = NULL;
		w->E = NULL;
	}
	if (initCone(k) < 0) {
		scs_printf("ERROR: initCone failure\n");
		scs_finish(d, w);
		return NULL;
	}
	w->p = initPriv(d);
	if (!w->p) {
		scs_printf("ERROR: initPriv failure\n");
		scs_finish(d, w);
		return NULL;
	}
	return w;
}

static void updateWork(Data * d, Work * w, Sol * sol) {
	/* before normalization */
	scs_int n = d->n;
	scs_int m = d->m;
	w->nm_b = calcNorm(d->b, m);
	w->nm_c = calcNorm(d->c, n);
#ifdef EXTRAVERBOSE
	printArray(d->b, d->m, "b");
	scs_printf("pre-normalized norm b = %4f\n", calcNorm(d->b, d->m));
	printArray(d->c, d->n, "c");
	scs_printf("pre-normalized norm c = %4f\n", calcNorm(d->c, d->n));
#endif
	if (d->normalize) {
		normalizeBC(d, w);
#ifdef EXTRAVERBOSE
		printArray(d->b, d->m, "bn");
		scs_printf("sc_b = %4f\n", w->sc_b);
		scs_printf("post-normalized norm b = %4f\n", calcNorm(d->b, d->m));
		printArray(d->c, d->n, "cn");
		scs_printf("sc_c = %4f\n", w->sc_c);
		scs_printf("post-normalized norm c = %4f\n", calcNorm(d->c, d->n));
#endif
	}
	if (d->warm_start) {
		warmStartVars(d, w, sol);
	} else {
		coldStartVars(d, w);
	}
	memcpy(w->h, d->c, n * sizeof(scs_float));
	memcpy(&(w->h[d->n]), d->b, m * sizeof(scs_float));
	memcpy(w->g, w->h, (n + m) * sizeof(scs_float));
	solveLinSys(d, w->p, w->g, NULL, -1);
	scaleArray(&(w->g[d->n]), -1, m);
	w->gTh = innerProd(w->h, w->g, n + m);
}

scs_int scs_solve(Work * w, Data * d, Cone * k, Sol * sol, Info * info) {
	scs_int i;
	timer solveTimer;
	struct residuals r;
	if (!d || !k || !sol || !info || !w || !d->b || !d->c) {
		scs_printf("ERROR: NULL input\n");
		return FAILURE;
	}
	tic(&solveTimer);
	info->statusVal = 0; /* not yet converged */
	updateWork(d, w, sol);
	if (d->verbose)
		printHeader(d, w, k);
	/* scs: */
	for (i = 0; i < d->max_iters; ++i) {
		memcpy(w->u_prev, w->u, (d->n + d->m + 1) * sizeof(scs_float));

		if (projectLinSys(d, w, i) < 0) return failureDefaultReturn(d, w, sol, info, "error in projectLinSys");
		if (projectCones(d, w, k, i) < 0) return failureDefaultReturn(d, w, sol, info, "error in projectCones");
		updateDualVars(d, w);

		if ((info->statusVal = converged(d, w, &r, i)) != 0)
			break;

		if (i % PRINT_INTERVAL == 0) {
			if (d->verbose) {
				printSummary(i, &r, &solveTimer);
#ifdef EXTRAVERBOSE
				 scs_printf("Norm u = %4f, ", calcNorm(w->u, d->n + d->m + 1));
				 scs_printf("Norm u_t = %4f, ", calcNorm(w->u_t, d->n + d->m + 1));
				 scs_printf("Norm v = %4f, ", calcNorm(w->v, d->n + d->m + 1));
				 scs_printf("tau = %4f, ", w->u[d->n + d->m]);
				 scs_printf("kappa = %4f, ", w->v[d->n + d->m]);
				 scs_printf("|u - u_prev| = %4f, ", calcNormDiff(w->u, w->u_prev, d->n + d->m + 1));
				 scs_printf("|u - u_t| = %4f\n", calcNormDiff(w->u, w->u_t, d->n + d->m + 1));
#endif
			}
		}
	}
	if (d->verbose) {
		printSummary(i, &r, &solveTimer);
	}
	setSolution(d, w, sol, info);
	/* populate info */
	info->iter = i;
	getInfo(d, w, sol, info);
	info->solveTime = tocq(&solveTimer);

	if (d->verbose)
		printFooter(d, w, info);
	/* un-normalize sol, b, c but not A */
	if (d->normalize)
		unNormalizeSolBC(d, w, sol);
	return info->statusVal;
}

void scs_finish(Data * d, Work * w) {
	finishCone();
	if (w) {
		if (d && d->normalize)
			unNormalizeA(d, w);
		freePriv(w->p);
		freeWork(w);
	}
}

Work * scs_init(Data * d, Cone * k, Info * info) {
	Work * w;
	timer initTimer;
	if (!d || !k || !info) {
		scs_printf("ERROR: Missing Data, Cone or Info input\n");
		return NULL;
	}
#ifdef EXTRAVERBOSE
	printData(d);
	printConeData(k);
#endif
#ifndef NOVALIDATE
	if (validate(d, k) < 0) {
		scs_printf("ERROR: Validation returned failure\n");
		return NULL;
	}
#endif
	tic(&initTimer);
	w = initWork(d, k);
	/* strtoc("init", &initTimer); */
	info->setupTime = tocq(&initTimer);
	if (d->verbose) {
		scs_printf("Setup time: %1.2es\n", info->setupTime / 1e3);
	}
	return w;
}

/* this just calls scs_init, scs_solve, and scs_finish */
scs_int scs(Data * d, Cone * k, Sol * sol, Info * info) {
#if ( defined _WIN32 || defined _WIN64 ) && !defined MATLAB_MEX_FILE && !defined PYTHON
	/* sets width of exponent for floating point numbers to 2 instead of 3 */
	unsigned int old_output_format = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
	Work * w = scs_init(d, k, info);
#ifdef EXTRAVERBOSE
	scs_printf("size of scs_int = %lu, size of scs_float = %lu\n", sizeof(scs_int), sizeof(scs_float));
#endif
	if (!w) {
		return failureDefaultReturn(d, NULL, sol, info, "could not initialize work");
	}
	scs_solve(w, d, k, sol, info);
	scs_finish(d, w);
	return info->statusVal;
}

/* TODO: re-integrate this eventually
 *  approximate convergence check:
 *  abs to prevent negative stopping tol */
/*
 scs_float tau = ABS(w->u[w->l-1]);
 scs_float kap = ABS(w->v[w->l-1]);
 r->resPri = calcNormDiff(w->u, w->u_t, w->l);
 r->resDual = calcNormDiff(w->u, w->u_prev, w->l);
 r->tau = tau;
 r->kap = kap;
 if (MIN(tau,kap)/MAX(tau,kap) < 1e-6 && MAX(r->resPri, r->resDual) < d->eps*(tau+kap)){
 return 1;
 }
 return 0
 */

