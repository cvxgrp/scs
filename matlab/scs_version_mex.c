#include "mex.h"
#include "scs.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* matlab usage: ver = scs_version() */
    if (nrhs != 0) {
        mexErrMsgTxt("scs_version doesn't take inout arguments.");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("scs_version returns 1 output argument only.");
    }
    plhs[0] = mxCreateString(scs_version());
    return;
}
