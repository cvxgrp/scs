#include "mex.h" 
#include "matrix.h"
#include "glbopts.h"
#include "scs.h"
#include "linAlg.h"

void freeMex(Data * d, Cone * k);

idxint parseWarmStart(const mxArray * p_mex, pfloat ** p, idxint l) {
    *p = scs_calloc(l, sizeof(pfloat));
    if (p_mex == NULL){
        return 0;
    } else if (mxIsSparse(p_mex) || (idxint)*mxGetDimensions(p_mex) != l) {
        scs_printf("Error parsing warm start input (make sure vectors are not sparse and of correct size), running without full warm-start");
        return 0;
    } else {
        memcpy(*p, mxGetPr(p_mex), l*sizeof(pfloat));
        return 1;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* matlab usage: scs(data,cone,params); */
  idxint i;  
  if (nrhs != 3){
    mexErrMsgTxt("Three arguments are required in this order: data struct, cone struct, params struct");
  }
  Data * d = mxMalloc(sizeof(Data)); 
  Cone * k = mxMalloc(sizeof(Cone)); 
  const mxArray *data = prhs[0];
   
  const mxArray *A_mex = (mxArray *) mxGetField(data,0,"A");
  if(A_mex == NULL) {
    scs_free(d); scs_free(k);
    mexErrMsgTxt("Data struct must contain a `A` entry.");
  }
  if (!mxIsSparse(A_mex)){
    scs_free(d); scs_free(k);
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  const mxArray *b_mex = (mxArray *) mxGetField(data,0,"b");
  if(b_mex == NULL) {
    scs_free(d); scs_free(k);
    mexErrMsgTxt("Data struct must contain a `b` entry.");
  }
  if (mxIsSparse(b_mex)){
    scs_free(d); scs_free(k);
    mexErrMsgTxt("Input vector b must be in dense format (pass in full(b))");
  }
  const mxArray *c_mex = (mxArray *) mxGetField(data,0,"c"); 
  if(c_mex == NULL) {
    scs_free(d); scs_free(k);
    mexErrMsgTxt("Data struct must contain a `c` entry.");
  }
  if (mxIsSparse(c_mex)){
    scs_free(d); scs_free(k);
    mexErrMsgTxt("Input vector c must be in dense format (pass in full(c))");
  }
  
  const mxArray *cone = prhs[1];
  const mxArray *params = prhs[2];
  d->n = (idxint)*(mxGetDimensions(c_mex));
  d->m = (idxint)*(mxGetDimensions(b_mex));
 
  d->b = mxGetPr(b_mex);
  d->c = mxGetPr(c_mex);
 
  const mxArray *ALPHA_mex = mxGetField(params,0,"ALPHA");
  if (ALPHA_mex == NULL) d->ALPHA = 1.8;
  else d->ALPHA = (pfloat)*mxGetPr(ALPHA_mex);
  
  const mxArray *RHO_X_mex = mxGetField(params,0,"RHO_X");
  if (RHO_X_mex == NULL) d->RHO_X = 1e-3;
  else d->RHO_X = (pfloat)*mxGetPr(RHO_X_mex);

  const mxArray *MAX_ITERS_mex = mxGetField(params,0,"MAX_ITERS");
  if (MAX_ITERS_mex == NULL) d->MAX_ITERS = 2500;
  else d->MAX_ITERS = (idxint)*mxGetPr(MAX_ITERS_mex);

  const mxArray *SCALE_mex = mxGetField(params,0,"SCALE");
  if (SCALE_mex == NULL) d->SCALE = 1;
  else d->SCALE = (pfloat)*mxGetPr(SCALE_mex);

  const mxArray *EPS_mex = mxGetField(params,0,"EPS");
  if (EPS_mex == NULL) d->EPS = 1e-3;
  else d->EPS = (pfloat)*mxGetPr(EPS_mex);

  const mxArray *UNDET_TOL_mex = mxGetField(params,0,"UNDET_TOL");
  if (UNDET_TOL_mex == NULL) d->UNDET_TOL = 1e-9;
  else d->UNDET_TOL = (pfloat)*mxGetPr(UNDET_TOL_mex);

  const mxArray *VERBOSE_mex = mxGetField(params,0,"VERBOSE");
  if (VERBOSE_mex == NULL) d->VERBOSE = 1;
  else d->VERBOSE = (idxint)*mxGetPr(VERBOSE_mex);

  const mxArray *NORMALIZE_mex = mxGetField(params,0,"NORMALIZE");
  if (NORMALIZE_mex == NULL) d->NORMALIZE = 1;
  else d->NORMALIZE = (idxint)*mxGetPr(NORMALIZE_mex);

  mxArray * kf = mxGetField(cone,0,"f");
  if(kf) k->f = (idxint)*mxGetPr(kf);
  else k->f = 0;
 
  mxArray * kl = mxGetField(cone,0,"l");
  if(kl) k->l = (idxint)*mxGetPr(kl);
  else k->l = 0;
  
  mxArray * ep = mxGetField(cone,0,"ep");
  if(ep) k->ep = (idxint)*mxGetPr(ep);
  else k->ep = 0;
 
  mxArray * ed = mxGetField(cone,0,"ed");
  if(ed) k->ed = (idxint)*mxGetPr(ed);
  else k->ed = 0;

  mxArray * kq = mxGetField(cone,0,"q");
  if (kq) {
      pfloat * q_mex = mxGetPr(kq);
      idxint ns = mxGetNumberOfDimensions(kq);
      const size_t * q_dims = mxGetDimensions(kq);
      k->qsize = q_dims[0];
      if (ns > 1 && q_dims[0] == 1) {
        k->qsize = q_dims[1];
      }
      k->q = mxMalloc(sizeof(idxint)*k->qsize);
      for ( i=0; i < k->qsize; i++ ){
          k->q[i] = (idxint)q_mex[i]; 
      }
  } else {
    k->qsize = 0;
    k->q = NULL;
  }

  mxArray * ks = mxGetField(cone,0,"s");
  if (ks) {
      pfloat * s_mex = mxGetPr(ks);
      idxint ns = mxGetNumberOfDimensions(ks);
      const size_t * s_dims = mxGetDimensions(ks);
      k->ssize = s_dims[0];
      if (ns > 1 && s_dims[0] == 1) {
        k->ssize = s_dims[1];
      }
      k->s = mxMalloc(sizeof(idxint)*k->ssize);
      for ( i=0; i < k->ssize; i++ ){
          k->s[i] = (idxint)s_mex[i]; 
      }
  } else {
      k->ssize = 0;
      k->s = NULL;
  }

  /*d->Anz = (idxint)mxGetNzmax(A_mex); */
  d->Ax = (pfloat *)mxGetPr(A_mex);
  /*
  d->Ap = (idxint *)mxMalloc(sizeof(int)*d->Anz);
  d->Ai = (idxint *)mxMalloc(sizeof(int)*d->Anz);
  mwSize * A_i = (mwSize*)mxGetIr(A_mex);
  for (i = 0; i < d->Anz; i++) {
    d->Ai[i] = (idxint)A_i[i];
  }
  mwSize * A_p = (mwSize*)mxGetJc(A_mex);
  for (i = 0; i < (d->n)+1; i++) {          
    d->Ap[i] = (idxint)A_p[i];
  }
  */

  /* XXX: 
  * these return (mwSize *), equivalent to (size_t *)
  * casting as (idxint *), when idxint = long seems to work
  * although maybe not on all machines:
  */
  d->Ap = (idxint *)mxGetJc(A_mex);
  d->Ai = (idxint *)mxGetIr(A_mex);
  
  /* warm-start inputs */
  Sol sol = { 0 };
  d->WARM_START = parseWarmStart((mxArray *) mxGetField(data,0,"x"), &(sol.x), d->n);
  d->WARM_START |= parseWarmStart((mxArray *) mxGetField(data,0,"y"), &(sol.y), d->m);
  d->WARM_START |= parseWarmStart((mxArray *) mxGetField(data,0,"s"), &(sol.s), d->m);
  
  Info info;
  idxint status = scs(d,k,&sol,&info);

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

  const char * infoFields[] = {"iter","status","pobj","dobj","resPri","resDual","relGap","setupTime","solveTime"}; 
  const idxint numInfoFields = 9;
  mwSize one[1] = {1};
  mxArray * xm;
  plhs[3] = mxCreateStructArray(1,one,numInfoFields,infoFields);

  mxSetField(plhs[3], 0, "status", mxCreateString(info.status));
   
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "iter", xm);
  *mxGetPr(xm) = info.iter;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "pobj", xm);
  *mxGetPr(xm) = info.pobj;

  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "dobj", xm);
  *mxGetPr(xm) = info.dobj;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "resPri", xm);
  *mxGetPr(xm) = info.resPri;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "resDual", xm);
  *mxGetPr(xm) = info.resDual;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "relGap", xm);
  *mxGetPr(xm) = info.relGap;

  /*info.time is millisecs - return value in secs */
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "setupTime", xm);
  *mxGetPr(xm) = info.setupTime/1e3; 

  /*info.time is millisecs - return value in secs */
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[3], 0, "solveTime", xm);
  *mxGetPr(xm) = info.solveTime/1e3; 

  freeMex(d, k);
  return; 
}

void freeMex(Data * d, Cone * k) {
  if(k->q) scs_free(k->q);
  if(k->s) scs_free(k->s);
  /* only needed if this data malloc-ed here */
  /*if(d->Ap) scs_free(d->Ap); */
  /*if(d->Ai) scs_free(d->Ai); */
  if(d) scs_free(d);
  if(k) scs_free(k);
}
