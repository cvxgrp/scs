#include "glbopts.h"
#include "scs.h"
#include "cones.h"
#include "linsys/amatrix.h"

#ifdef INDIRECTJ
#include "org_scs_IndirectSolver.h"
#else
#include "org_scs_DirectSolver.h"
#endif

jobject getObjUsingGetter(JNIEnv * env, jobject obj, char * method, char * signature) {
    // Get method ID for method getSomeArray that returns an array
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, signature);
    // Call the method, returns JObject (because Array is instance of Object)
    return (*env)->CallObjectMethod(env, obj, mid);
}

scs_int getIntUsingGetter(JNIEnv * env, jobject obj, char * method) {
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "()I");
    return (scs_int) (*env)->CallIntMethod(env, obj, mid);
}

scs_float getFloatUsingGetter(JNIEnv * env, jobject obj, char * method) {
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "()D");
    return (scs_float) (*env)->CallDoubleMethod(env, obj, mid);
}

scs_int getBooleanUsingGetter(JNIEnv * env, jobject obj, char * method) {
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "()Z");
    return (scs_int) (*env)->CallBooleanMethod(env, obj, mid);
}

scs_int * getIntArrayUsingGetter(JNIEnv * env, jobject obj, char * method, scs_int * length) {
    scs_int * out;
    jintArray arr = (jintArray) getObjUsingGetter(env, obj, method, "()[I");
    if (!arr) {
        *length = 0;
        return SCS_NULL;
    }
    *length = (*env)->GetArrayLength(env, arr);
    out = scs_malloc(sizeof(scs_int) * (*length));
    (*env)->GetIntArrayRegion( env, arr, 0, *length, out);
    return out;
}

scs_float * getFloatArrayUsingGetter(JNIEnv * env, jobject obj, char * method, scs_int * length) {
    scs_float * out;
    jdoubleArray arr = (jdoubleArray) getObjUsingGetter(env, obj, method, "()[D");
    if (!arr) {
        *length = 0;
        return SCS_NULL;
    }
    *length = (*env)->GetArrayLength(env, arr);
    out = scs_malloc(sizeof(scs_float) * (*length));
    (*env)->GetDoubleArrayRegion( env, arr, 0, *length, out);
    return out;
}

Cone * getConeStruct(JNIEnv * env, jobject coneJava) {
    Cone * k = scs_calloc(1, sizeof(Cone));
    k->q = getIntArrayUsingGetter(env, coneJava, "getQ", &(k->qsize));
    k->s = getIntArrayUsingGetter(env, coneJava, "getS", &(k->ssize));
    k->l = getIntUsingGetter(env, coneJava, "getL");
    k->f = getIntUsingGetter(env, coneJava, "getF");
    k->ep = getIntUsingGetter(env, coneJava, "getEp");
    k->ed = getIntUsingGetter(env, coneJava, "getEd");
    k->p = getFloatArrayUsingGetter(env, coneJava, "getP", &(k->psize));
    return k;
}

AMatrix * getAMatrix(JNIEnv *env, jobject AJava, scs_int m, scs_int n) {
    scs_int leni, lenp, lenx;
    AMatrix * A = scs_calloc(1, sizeof(AMatrix));
    // populate A
    A->i = getIntArrayUsingGetter(env, AJava, "getRowIdxs", &leni);
    A->p = getIntArrayUsingGetter(env, AJava, "getColIdxs", &lenp);
    A->x = getFloatArrayUsingGetter(env, AJava, "getValues", &lenx);
    A->m = m;
    A->n = n;
    return A;
}

void populateParams(JNIEnv * env, jobject paramsJava, Data * d) {
    d->stgs = scs_malloc(sizeof(Settings));
    d->stgs->max_iters = getIntUsingGetter(env, paramsJava, "getMaxIters");
    d->stgs->eps = getFloatUsingGetter(env, paramsJava, "getEps");
    d->stgs->alpha = getFloatUsingGetter(env, paramsJava, "getAlpha");
    d->stgs->rho_x = getFloatUsingGetter(env, paramsJava, "getRhoX");
    d->stgs->cg_rate = getFloatUsingGetter(env, paramsJava, "getCgRate");
    d->stgs->verbose = getBooleanUsingGetter(env, paramsJava, "isVerbose");
    d->stgs->normalize = getBooleanUsingGetter(env, paramsJava, "isNormalize");
    d->stgs->scale = getFloatUsingGetter(env, paramsJava, "getScale");
    d->stgs->warm_start = getBooleanUsingGetter(env, paramsJava, "isWarmStart");
}

Data * getDataStruct(JNIEnv * env, jobject AJava, jdoubleArray bJava, jdoubleArray cJava, jobject paramsJava) {
    Data * d = scs_calloc(1, sizeof(Data));
    d->b = (*env)->GetDoubleArrayElements(env, bJava, SCS_NULL);
    d->m = (*env)->GetArrayLength(env, bJava);
    d->c = (*env)->GetDoubleArrayElements(env, cJava, SCS_NULL);
    d->n = (*env)->GetArrayLength(env, cJava);
    d->A = getAMatrix(env, AJava, d->m, d->n);
    populateParams(env, paramsJava, d);
    return d;
}

void setFloatArrayUsingSetter(JNIEnv * env, jobject obj, scs_float * arr, scs_int length, char * method) {
    jdoubleArray out = (*env)->NewDoubleArray(env, length);
    (*env)->SetDoubleArrayRegion(env, out, 0, length, arr);
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "([D)V");
    (*env)->CallVoidMethod(env, obj, mid, out);
}

void setStringUsingSetter(JNIEnv * env, jobject obj, char * str, char * method) {
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "(Ljava/lang/String;)V");
    (*env)->CallVoidMethod(env, obj, mid, (*env)->NewStringUTF(env, str));
}

void setIntUsingSetter(JNIEnv * env, jobject obj, scs_int i, char * method) {
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "(I)V");
    (*env)->CallVoidMethod(env, obj, mid, i);
}

void setFloatUsingSetter(JNIEnv * env, jobject obj, scs_float f, char * method) {
    jclass clazz = (*env)->GetObjectClass(env, obj);
    jmethodID mid = (*env)->GetMethodID(env, clazz, method, "(D)V");
    (*env)->CallVoidMethod(env, obj, mid, f);
}

void setSol(JNIEnv * env, jobject solJava, Data * d, Sol * sol) {
    setFloatArrayUsingSetter(env, solJava, sol->x, d->n, "setX");
    setFloatArrayUsingSetter(env, solJava, sol->y, d->m, "setY");
    setFloatArrayUsingSetter(env, solJava, sol->s, d->m, "setS");
}

void setInfo(JNIEnv * env, jobject infoJava, Info * info) {
    setIntUsingSetter(env, infoJava, info->iter, "setIter");
    setIntUsingSetter(env, infoJava, info->statusVal, "setStatusVal");
    setStringUsingSetter(env, infoJava, info->status, "setStatus");
    setFloatUsingSetter(env, infoJava, info->pobj, "setPobj");
    setFloatUsingSetter(env, infoJava, info->dobj, "setDobj");
    setFloatUsingSetter(env, infoJava, info->resPri, "setResPri");
    setFloatUsingSetter(env, infoJava, info->resDual, "setResDual");
    setFloatUsingSetter(env, infoJava, info->resInfeas, "setResInfeas");
    setFloatUsingSetter(env, infoJava, info->resUnbdd, "setResUnbdd");
    setFloatUsingSetter(env, infoJava, info->relGap, "setRelGap");
    setFloatUsingSetter(env, infoJava, info->setupTime, "setSetupTime");
    setFloatUsingSetter(env, infoJava, info->solveTime, "setSolveTime");
}

#ifdef INDIRECTJ
JNIEXPORT jstring JNICALL Java_org_scs_IndirectSolver_cversion(JNIEnv *env, jclass clazz)
#else
JNIEXPORT jstring JNICALL Java_org_scs_DirectSolver_cversion(JNIEnv *env, jclass clazz)
#endif
{
    return (*env)->NewStringUTF(env, scs_version());
}

#ifdef INDIRECTJ
JNIEXPORT void JNICALL Java_org_scs_IndirectSolver_csolve(JNIEnv *env, jclass clazz, jobject AJava,
        jdoubleArray bJava, jdoubleArray cJava, jobject coneJava, jobject paramsJava, jobject solJava, jobject infoJava)
#else
JNIEXPORT void JNICALL Java_org_scs_DirectSolver_csolve(JNIEnv *env, jclass clazz, jobject AJava,
        jdoubleArray bJava, jdoubleArray cJava, jobject coneJava, jobject paramsJava, jobject solJava, jobject infoJava)
#endif
{
    /* Parse out the data into C form, then pass to SCS, the convert solution back to java object */
    /* Assume AJava contains matrix in column compressed sparse format */
    Data * d = getDataStruct(env, AJava, bJava, cJava, paramsJava);
    Cone * k = getConeStruct(env, coneJava);

    Sol * sol = scs_calloc(1, sizeof(Sol));
    Info * info = scs_calloc(1, sizeof(Info));
    scs(d, k, sol, info);

    setSol(env, solJava, d, sol);
    setInfo(env, infoJava, info);

    freeData(d, k);
    freeSol(sol);
    scs_free(info);
}

