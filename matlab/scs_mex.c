#include "mex.h"
#include "matrix.h"
#include "glbopts.h"
#include "scs.h"
#include "linAlg.h"
#include "linsys/amatrix.h"

void freeMex(Data * d, Cone * k);

idxint parseWarmStart(const mxArray * p_mex, pfloat ** p, idxint l) {
	*p = scs_calloc(l, sizeof(pfloat)); /* this allocates memory used for Sol */
	if (p_mex == NULL) {
		return 0;
	} else if (mxIsSparse(p_mex) || (idxint) *mxGetDimensions(p_mex) != l) {
		scs_printf("Error parsing warm start input (make sure vectors are not sparse and of correct size), running without full warm-start");
		return 0;
	} else {
		memcpy(*p, mxGetPr(p_mex), l * sizeof(pfloat));
		return 1;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* matlab usage: scs(data,cone,params); */
	idxint i, ns, status;
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
	const pfloat *q_mex;
	const pfloat *s_mex;
	const size_t *q_dims;
	const size_t *s_dims;

	const mxArray *cone;
	const mxArray *params;

	const mwSize one[1] = { 1 };
	const int numInfoFields = 9;
	const char * infoFields[] = { "iter", "status", "pobj", "dobj", "resPri", "resDual", "relGap", "setupTime",
			"solveTime" };
	mxArray *tmp;


	if (nrhs != 3) {
		mexErrMsgTxt("Three arguments are required in this order: data struct, cone struct, params struct");
	}
	if (nlhs > 4) {
		mexErrMsgTxt("scs returns up to 4 output arguments only.");
	}
	d = mxMalloc(sizeof(Data));
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
	params = prhs[2];
	d->n = (idxint) *(mxGetDimensions(c_mex));
	d->m = (idxint) *(mxGetDimensions(b_mex));

	d->b = (pfloat *)mxGetPr(b_mex);
	d->c = (pfloat *)mxGetPr(c_mex);

	/* params */
	tmp = mxGetField(params, 0, "ALPHA");
	if (tmp == NULL)
		d->ALPHA = 1.8;
	else
		d->ALPHA = (pfloat) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "RHO_X");
	if (tmp == NULL)
		d->RHO_X = 1e-3;
	else
		d->RHO_X = (pfloat) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "MAX_ITERS");
	if (tmp == NULL)
		d->MAX_ITERS = 2500;
	else
		d->MAX_ITERS = (idxint) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "SCALE");
	if (tmp == NULL)
		d->SCALE = 5;
	else
		d->SCALE = (pfloat) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "EPS");
	if (tmp == NULL)
		d->EPS = 1e-3;
	else
		d->EPS = (pfloat) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "CG_RATE");
	if (tmp == NULL)
		d->CG_RATE = 2;
	else
		d->CG_RATE = (pfloat) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "VERBOSE");
	if (tmp == NULL)
		d->VERBOSE = 1;
	else
		d->VERBOSE = (idxint) *mxGetPr(tmp);

	tmp = mxGetField(params, 0, "NORMALIZE");
	if (tmp == NULL)
		d->NORMALIZE = 1;
	else
		d->NORMALIZE = (idxint) *mxGetPr(tmp);

	/* cones */
	kf = mxGetField(cone, 0, "f");
	if (kf && !mxIsEmpty(kf))
		k->f = (idxint) *mxGetPr(kf);
	else
		k->f = 0;

	kl = mxGetField(cone, 0, "l");
	if (kl && !mxIsEmpty(kl))
		k->l = (idxint) *mxGetPr(kl);
	else
		k->l = 0;

	kep = mxGetField(cone, 0, "ep");
	if (kep && !mxIsEmpty(kep))
		k->ep = (idxint) *mxGetPr(kep);
	else
		k->ep = 0;

	ked = mxGetField(cone, 0, "ed");
	if (ked && !mxIsEmpty(ked))
		k->ed = (idxint) *mxGetPr(ked);
	else
		k->ed = 0;

	kq = mxGetField(cone, 0, "q");
	if (kq && !mxIsEmpty(kq)) {
		q_mex = mxGetPr(kq);
		ns = (idxint) mxGetNumberOfDimensions(kq);
		q_dims = mxGetDimensions(kq);
		k->qsize = (idxint) q_dims[0];
		if (ns > 1 && q_dims[0] == 1) {
			k->qsize = (idxint) q_dims[1];
		}
		k->q = mxMalloc(sizeof(idxint) * k->qsize);
		for (i = 0; i < k->qsize; i++) {
			k->q[i] = (idxint) q_mex[i];
		}
	} else {
		k->qsize = 0;
		k->q = NULL;
	}

	ks = mxGetField(cone, 0, "s");
	if (ks && !mxIsEmpty(ks)) {
		s_mex = mxGetPr(ks);
		ns = (idxint) mxGetNumberOfDimensions(ks);
		s_dims = mxGetDimensions(ks);
		k->ssize = (idxint) s_dims[0];
		if (ns > 1 && s_dims[0] == 1) {
			k->ssize = (idxint) s_dims[1];
		}
		k->s = mxMalloc(sizeof(idxint) * k->ssize);
		for (i = 0; i < k->ssize; i++) {
			k->s[i] = (idxint) s_mex[i];
		}
	} else {
		k->ssize = 0;
		k->s = NULL;
	}
	A = scs_malloc(sizeof(AMatrix));
	A->x = (pfloat *) mxGetPr(A_mex);
	/* XXX:
	 * these return (mwIndex *), equivalent to (size_t *)
	 * casting as (idxint *), when idxint = long seems to work
	 * although maybe not on all machines:
	 */
	A->p = (idxint *) mxGetJc(A_mex);
	A->i = (idxint *) mxGetIr(A_mex);
	d->A = A;
	/* warm-start inputs, allocates sol->x, ->y, ->s even if warm start not used */
	d->WARM_START = parseWarmStart((mxArray *) mxGetField(data, 0, "x"), &(sol.x), d->n);
	d->WARM_START |= parseWarmStart((mxArray *) mxGetField(data, 0, "y"), &(sol.y), d->m);
	d->WARM_START |= parseWarmStart((mxArray *) mxGetField(data, 0, "s"), &(sol.s), d->m);

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
	*mxGetPr(tmp) = (pfloat) info.iter;

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
	mxSetField(plhs[3], 0, "relGap", tmp);
	*mxGetPr(tmp) = info.relGap;

	/*info.time is millisecs - return value in secs */
	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "setupTime", tmp);
	*mxGetPr(tmp) = info.setupTime / 1e3;

	/*info.time is millisecs - return value in secs */
	tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[3], 0, "solveTime", tmp);
	*mxGetPr(tmp) = info.solveTime / 1e3;

	freeMex(d, k);
	return;
}

void freeMex(Data * d, Cone * k) {
	if (k->q)
		scs_free(k->q);
	if (k->s)
		scs_free(k->s);
	if (d) {
		if(d->A) scs_free(d->A);
		scs_free(d);
	}
	if (k)
		scs_free(k);
}
