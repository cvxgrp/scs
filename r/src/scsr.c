#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "../../include/scs.h"
#include "../../include/util.h"
#include "../../linsys/amatrix.h"

SEXP getListElement(SEXP list, const char *str) {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < length(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    }
    return elmt;
}

SEXP floatVec2R(scs_int n, scs_float *in) {
    SEXP ret;
    scs_int i;
    scs_float *vec;

    PROTECT(ret = allocVector(REALSXP, n));
    vec = REAL(ret);
    for (i = 0; i < n; i++) {
        vec[i] = in[i];
    }
    return ret;
}

scs_float getFloatFromListWithDefault(SEXP list, const char *str,
                                      scs_float def) {
    SEXP val = getListElement(list, str);
    if (val == R_NilValue) {
        return def;
    }
    val = coerceVector(val, REALSXP);
    return REAL(val)[0];
}

scs_int getIntFromListWithDefault(SEXP list, const char *str, scs_int def) {
    SEXP val = getListElement(list, str);
    if (val == R_NilValue) {
        return def;
    }
    val = coerceVector(val, INTSXP);
    return INTEGER(val)[0];
}

scs_float *getFloatVectorFromList(SEXP list, const char *str, scs_int *len) {
    SEXP vec = getListElement(list, str);
    *len = length(vec);
    vec = coerceVector(vec, REALSXP);
    return REAL(vec);
}

scs_int *getIntVectorFromList(SEXP list, const char *str, scs_int *len) {
    SEXP vec = getListElement(list, str);
    *len = length(vec);
    vec = coerceVector(vec, INTSXP);
    return INTEGER(vec);
}

SEXP populateInfoR(Info *info, scs_int *num_protectaddr) {
    scs_int num_protected = 0;
    SEXP infor, info_names, iter_r, status_r, statusVal_r, pobj_r, dobj_r,
        resPri_r, resDual_r, resInfeas_r, resUnbdd_r, relGap_r, setupTime_r,
        solveTime_r;

    PROTECT(infor = NEW_LIST(12));
    num_protected++;
    PROTECT(info_names = NEW_CHARACTER(12));
    num_protected++;
    SET_NAMES(infor, info_names);

    PROTECT(iter_r = allocVector(INTSXP, 1));
    num_protected++;
    INTEGER(iter_r)[0] = info->iter;
    SET_STRING_ELT(info_names, 0, mkChar("iter"));
    SET_VECTOR_ELT(infor, 0, iter_r);

    PROTECT(status_r = NEW_CHARACTER(1));
    num_protected++;
    SET_STRING_ELT(status_r, 0, mkChar(info->status));
    SET_STRING_ELT(info_names, 1, mkChar("status"));
    SET_VECTOR_ELT(infor, 1, status_r);

    PROTECT(statusVal_r = allocVector(INTSXP, 1));
    num_protected++;
    INTEGER(statusVal_r)[0] = info->statusVal;
    SET_STRING_ELT(info_names, 2, mkChar("statusVal"));
    SET_VECTOR_ELT(infor, 2, statusVal_r);

    PROTECT(pobj_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(pobj_r)[0] = info->pobj;
    SET_STRING_ELT(info_names, 3, mkChar("pobj"));
    SET_VECTOR_ELT(infor, 3, pobj_r);

    PROTECT(dobj_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(dobj_r)[0] = info->dobj;
    SET_STRING_ELT(info_names, 4, mkChar("dobj"));
    SET_VECTOR_ELT(infor, 4, dobj_r);

    PROTECT(resPri_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(resPri_r)[0] = info->resPri;
    SET_STRING_ELT(info_names, 5, mkChar("resPri"));
    SET_VECTOR_ELT(infor, 5, resPri_r);

    PROTECT(resDual_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(resDual_r)[0] = info->resDual;
    SET_STRING_ELT(info_names, 6, mkChar("resDual"));
    SET_VECTOR_ELT(infor, 6, resDual_r);

    PROTECT(resInfeas_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(resInfeas_r)[0] = info->resInfeas;
    SET_STRING_ELT(info_names, 7, mkChar("resInfeas"));
    SET_VECTOR_ELT(infor, 7, resInfeas_r);

    PROTECT(resUnbdd_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(resUnbdd_r)[0] = info->resUnbdd;
    SET_STRING_ELT(info_names, 8, mkChar("resUnbdd"));
    SET_VECTOR_ELT(infor, 8, resUnbdd_r);

    PROTECT(relGap_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(relGap_r)[0] = info->relGap;
    SET_STRING_ELT(info_names, 9, mkChar("relGap"));
    SET_VECTOR_ELT(infor, 9, relGap_r);

    PROTECT(setupTime_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(setupTime_r)[0] = info->setupTime;
    SET_STRING_ELT(info_names, 10, mkChar("setupTime"));
    SET_VECTOR_ELT(infor, 10, setupTime_r);

    PROTECT(solveTime_r = allocVector(REALSXP, 1));
    num_protected++;
    REAL(solveTime_r)[0] = info->solveTime;
    SET_STRING_ELT(info_names, 11, mkChar("solveTime"));
    SET_VECTOR_ELT(infor, 11, solveTime_r);

    *num_protectaddr += num_protected;
    return infor;
}

SEXP scsr(SEXP data, SEXP cone, SEXP params) {
    scs_int len, num_protected = 0;
    SEXP ret, retnames, infor, xr, yr, sr;

    /* allocate memory */
    Data *d = scs_malloc(sizeof(Data));
    Cone *k = scs_malloc(sizeof(Cone));
    Settings *stgs = scs_malloc(sizeof(Settings));
    AMatrix *A = scs_malloc(sizeof(AMatrix));
    Info *info = scs_calloc(1, sizeof(Info));
    Sol *sol = scs_calloc(1, sizeof(Sol));

    d->b = getFloatVectorFromList(data, "b", &len);
    d->c = getFloatVectorFromList(data, "c", &len);
    d->n = getIntFromListWithDefault(data, "n", 0);
    d->m = getIntFromListWithDefault(data, "m", 0);

    A->m = d->m;
    A->n = d->n;
    A->x = getFloatVectorFromList(data, "Ax", &len);
    A->i = getIntVectorFromList(data, "Ai", &len);
    A->p = getIntVectorFromList(data, "Ap", &len);
    d->A = A;

    stgs->max_iters = getIntFromListWithDefault(params, "max_iters", MAX_ITERS);
    stgs->normalize = getIntFromListWithDefault(params, "normalize", NORMALIZE);
    stgs->verbose = getIntFromListWithDefault(params, "verbose", VERBOSE);
    stgs->cg_rate = getFloatFromListWithDefault(params, "cg_rate", CG_RATE);
    stgs->scale = getFloatFromListWithDefault(params, "scale", SCALE);
    stgs->rho_x = getFloatFromListWithDefault(params, "rho_x", RHO_X);
    stgs->alpha = getFloatFromListWithDefault(params, "alpha", ALPHA);
    stgs->eps = getFloatFromListWithDefault(params, "eps", EPS);
    /* TODO add warm starting */
    stgs->warm_start =
        getIntFromListWithDefault(params, "warm_start", WARM_START);
    d->stgs = stgs;

    k->f = getIntFromListWithDefault(cone, "f", 0);
    k->l = getIntFromListWithDefault(cone, "l", 0);
    k->ep = getIntFromListWithDefault(cone, "ep", 0);
    k->ed = getIntFromListWithDefault(cone, "ed", 0);
    k->q = getIntVectorFromList(cone, "q", &(k->qsize));
    k->s = getIntVectorFromList(cone, "s", &(k->ssize));
    k->p = getFloatVectorFromList(cone, "p", &(k->psize));

    /* solve! */
    scs(d, k, sol, info);

    infor = populateInfoR(info, &num_protected);

    PROTECT(ret = NEW_LIST(4));
    num_protected++;
    PROTECT(retnames = NEW_CHARACTER(4));
    num_protected++;
    SET_NAMES(ret, retnames);

    xr = floatVec2R(d->n, sol->x);
    num_protected++;
    yr = floatVec2R(d->m, sol->y);
    num_protected++;
    sr = floatVec2R(d->m, sol->s);
    num_protected++;

    SET_STRING_ELT(retnames, 0, mkChar("x"));
    SET_VECTOR_ELT(ret, 0, xr);
    SET_STRING_ELT(retnames, 1, mkChar("y"));
    SET_VECTOR_ELT(ret, 1, yr);
    SET_STRING_ELT(retnames, 2, mkChar("s"));
    SET_VECTOR_ELT(ret, 2, sr);
    SET_STRING_ELT(retnames, 3, mkChar("info"));
    SET_VECTOR_ELT(ret, 3, infor);

    /* free memory */
    scs_free(info);
    scs_free(d);
    scs_free(k);
    scs_free(stgs);
    scs_free(A);
    freeSol(sol);
    UNPROTECT(num_protected);
    return ret;
}
