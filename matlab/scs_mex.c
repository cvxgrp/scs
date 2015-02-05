#include "mex.h"
#include "matrix.h"
#include "glbopts.h"
#include "scs.h"
#include "linAlg.h"
#include "linsys/amatrix.h"

void freeMex(Data * d, Cone * k);

scs_int parseWarmStart(const mxArray * p_mex, scs_float ** p, scs_int l) {
	*p = scs_calloc(l, sizeof(scs_float)); /* this allocates memory used for Sol */
	if (p_mex == NULL) {
		return 0;
	} else if (mxIsSparse(p_mex) || (scs_int) *mxGetDimensions(p_mex) != l) {
		scs_printf("Error parsing warm start input (make sure vectors are not sparse and of correct size), running without full warm-start");
		return 0;
	} else {
		memcpy(*p, mxGetPr(p_mex), l * sizeof(scs_float));
		return 1;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* matlab usage: scs(data,cone,settings); */
	scs_int i, ns, status;
	Data *d;
	Cone *k;
	Sol sol = { 0 };
	Info info;
	AMatrix * A;

	const mxArray *data;
	const mxArray *A_mex;
	const mxArray *b_mex;
	const mxArray *c_mex;

	const mxArray *kf;
	const mxArray *kl;
	const mxArray *kq;
	const mxArray *ks;
	const mxArray *kep;
	const mxArray *ked;
	const mxArray *kp;
	const scs_float *q_mex;
    const scs_float *s_mex;
    const scs_float *p_mex;
    const size_t *q_dims;
    const size_t *s_dims;
    const size_t *p_dims;

	const mxArray *cone;
	const mxArray *settings;

	const mwSize one[1] = { 1 };
	const int numInfoFields = 11;
	const char * infoFields[] = { "iter", "status", "pobj", "dobj", "resPri", "resDual", "resInfeas", "resUnbdd",
		"relGap", "setupTime", "solveTime" };
	mxArray *tmp;


	if (nrhs != 3) {
		mexErrMsgTxt("Three arguments are required in this order: data struct, cone struct, settings struct");
	}
	if (nlhs > 4) {
		mexErrMsgTxt("scs returns up to 4 output arguments only.");
	}
	d = mxMalloc(sizeof(Data));
    d->stgs = mxMalloc(sizeof(Settings));
	k = mxMalloc(sizeof(Cone));
	data = prhs[0];

	A_mex = (mxArray *) mxGetField(data, 0, "A");
	if (A_mex == NULL) {
		scs_free(d);
		scs_free(k);
		mexErrMsgTxt("Data struct must contain a `A` entry.");
	}
	if (!mxIsSparse(A_mex)) {
		scs_free(d);
		scs_free(k);
		mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
	}
	b_mex = (mxArray *) mxGetField(data, 0, "b");
	if (b_mex == NULL) {
		scs_free(d);
		scs_free(k);
		mexErrMsgTxt("Data struct must contain a `b` entry.");
	}
	if (mxIsSparse(b_mex)) {
		scs_free(d);
		scs_free(k);
		mexErrMsgTxt("Input vector b must be in dense format (pass in full(b))");
	}
	c_mex = (mxArray *) mxGetField(data, 0, "c");
	if (c_mex == NULL) {
		scs_free(d);
		scs_free(k);
		mexErrMsgTxt("Data struct must contain a `c` entry.");
	}
	if (mxIsSparse(c_mex)) {
		scs_free(d);
		scs_free(k);
		mexErrMsgTxt("Input vector c must be in dense format (pass in full(c))");
	}

	cone = prhs[1];
	settings = prhs[2];
	d->n = (scs_int) *(mxGetDimensions(c_mex));
	d->m = (scs_int) *(mxGetDimensions(b_mex));

	d->b = (scs_float *)mxGetPr(b_mex);
	d->c = (scs_float *)mxGetPr(c_mex);
	setDefaultSettings(d);

	/* settings */
	tmp = mxGetField(settings, 0, "alpha");
	if (tmp != NULL)
		d->stgs->alpha = (scs_float) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "rho_x");
	if (tmp != NULL)
		d->stgs->rho_x = (scs_float) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "max_iters");
	if (tmp != NULL)
		d->stgs->max_iters = (scs_int) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "scale");
	if (tmp != NULL)
		d->stgs->scale = (scs_float) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "eps");
	if (tmp != NULL)
		d->stgs->eps = (scs_float) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "cg_rate");
	if (tmp != NULL)
		d->stgs->cg_rate = (scs_float) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "verbose");
	if (tmp != NULL)
		d->stgs->verbose = (scs_int) *mxGetPr(tmp);

	tmp = mxGetField(settings, 0, "normalize");
	if (tmp != NULL)
		d->stgs->normalize = (scs_int) *mxGetPr(tmp);

	/* cones */
	kf = mxGetField(cone, 0, "f");
	if (kf && !mxIsEmpty(kf))
		k->f = (scs_int) *mxGetPr(kf);
	else
		k->f = 0;

	kl = mxGetField(cone, 0, "l");
	if (kl && !mxIsEmpty(kl))
		k->l = (scs_int) *mxGetPr(kl);
	else
		k->l = 0;

	kep = mxGetField(cone, 0, "ep");
	if (kep && !mxIsEmpty(kep))
		k->ep = (scs_int) *mxGetPr(kep);
	else
		k->ep = 0;

	ked = mxGetField(cone, 0, "ed");
	if (ked && !mxIsEmpty(ked))
		k->ed = (scs_int) *mxGetPr(ked);
	else
		k->ed = 0;

	kq = mxGetField(cone, 0, "q");
	if (kq && !mxIsEmpty(kq)) {
		q_mex = mxGetPr(kq);
		ns = (scs_int) mxGetNumberOfDimensions(kq);
		q_dims = mxGetDimensions(kq);
		k->qsize = (scs_int) q_dims[0];
		if (ns > 1 && q_dims[0] == 1) {
			k->qsize = (scs_int) q_dims[1];
		}
		k->q = mxMalloc(sizeof(scs_int) * k->qsize);
		for (i = 0; i < k->qsize; i++) {
			k->q[i] = (scs_int) q_mex[i];
		}
	} else {
		k->qsize = 0;
		k->q = NULL;
	}

	ks = mxGetField(cone, 0, "s");
	if (ks && !mxIsEmpty(ks)) {
		s_mex = mxGetPr(ks);
		ns = (scs_int) mxGetNumberOfDimensions(ks);
		s_dims = mxGetDimensions(ks);
		k->ssize = (scs_int) s_dims[0];
		if (ns > 1 && s_dims[0] == 1) {
			k->ssize = (scs_int) s_dims[1];
		}
		k->s = mxMalloc(sizeof(scs_int) * k->ssize);
		for (i = 0; i < k->ssize; i++) {
			k->s[i] = (scs_int) s_mex[i];
		}
	} else {
		k->ssize = 0;
		k->s = NULL;
	}

    kp = mxGetField(cone, 0, "p");
    if (kp && !mxIsEmpty(kp)) {
        p_mex = mxGetPr(kp);
        ns = (scs_int) mxGetNumberOfDimensions(kp);
        p_dims = mxGetDimensions(kp);
        k->psize = (scs_int) p_dims[0];
        if (ns > 1 && p_dims[0] == 1) {
            k->psize = (scs_int) p_dims[1];
        }
        k->p = mxMalloc(sizeof(scs_float) * k->psize);
        for (i = 0; i < k->psize; i++) {
            k->p[i] = (scs_float) p_mex[i];
        }
    } else {
        k->psize = 0;
        k->p = NULL;
    }

    A = scs_malloc(sizeof(AMatrix));
    A->n = d->n;
    A->m = d->m;
    A->x = (scs_float *) mxGetPr(A_mex);
    /* TODO:
     * these return (mwIndex *), equivalent to (size_t *)
     * casting as (scs_int *), when scs_int = long seems to work
     * although maybe not on all machines:
     */
    A->p = (scs_int *) mxGetJc(A_mex);
    A->i = (scs_int *) mxGetIr(A_mex);
    d->A = A;
    /* warm-start inputs, allocates sol->x, ->y, ->s even if warm start not used */
    d->stgs->warm_start = parseWarmStart((mxArray *) mxGetField(data, 0, "x"), &(sol.x), d->n);
    d->stgs->warm_start |= parseWarmStart((mxArray *) mxGetField(data, 0, "y"), &(sol.y), d->m);
    d->stgs->warm_start |= parseWarmStart((mxArray *) mxGetField(data, 0, "s"), &(sol.s), d->m);

	status = scs(d, k, &sol, &info);

	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetPr(plhs[0], sol.x);
	mxSetM(plhs[0], d->n);
	mxSetN(plhs[0], 1);

	plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetPr(plhs[1], sol.y);
	mxSetM(plhs[1], d->m);
	mxSetN(plhs[1], 1);

	plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetPr(plhs[2], sol.s);
	mxSetM(plhs[2], d->m);
	mxSetN(plhs[2], 1);

	plhs[3] = mxCreateStructArray(1, one, numInfoFields, infoFields);

	mxSetField(plhs[3], 0, "status", mxCreateString(info.status));

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "iter", tmp);
	*mxGetPr(tmp) = (scs_float) info.iter;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "pobj", tmp);
	*mxGetPr(tmp) = info.pobj;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "dobj", tmp);
	*mxGetPr(tmp) = info.dobj;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "resPri", tmp);
	*mxGetPr(tmp) = info.resPri;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "resDual", tmp);
	*mxGetPr(tmp) = info.resDual;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "resInfeas", tmp);
	*mxGetPr(tmp) = info.resInfeas;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "resUnbdd", tmp);
	*mxGetPr(tmp) = info.resUnbdd;

	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "relGap", tmp);
	*mxGetPr(tmp) = info.relGap;

	/*info.time is millisecs - return value in secs */
	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "setupTime", tmp);
	*mxGetPr(tmp) = info.setupTime;

	/*info.time is millisecs - return value in secs */
	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "solveTime", tmp);
	*mxGetPr(tmp) = info.solveTime;

	freeMex(d, k);
	return;
}

void freeMex(Data * d, Cone * k) {
	if (k->q)
		scs_free(k->q);
	if (k->s)
        scs_free(k->s);
    if (k->p)
        scs_free(k->p);
    if (d) {
        if(d->A) scs_free(d->A);
        if(d->stgs) scs_free(d->stgs);
        scs_free(d);
    }
    if (k)
        scs_free(k);
}
