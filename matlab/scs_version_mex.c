#include "mex.h"
#include "scs.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* matlab usage: ver = scs_version() */
    if (nrhs != 0) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    plhs[0] = mxCreateString(scs_version());
    return;
}
