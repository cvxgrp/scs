#include "glbopts.h"
#include "scs.h"
#include "cones.h"
#include "linsys/amatrix.h"

#ifdef INDIRECTJ
#include "scs_IndirectSolver.h"
#else
#include "scs_DirectSolver.h"
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
        return NULL;
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
        return NULL;
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
    return k;
}

AMatrix * getAMatrix(JNIEnv *env, jobject AJava) {
    scs_int leni, lenp, lenx;
    AMatrix * A = scs_calloc(1, sizeof(AMatrix));
    // populate A
    A->x = getFloatArrayUsingGetter(env, AJava, "getValues", &lenx);
    return A;
}

void populateParams(JNIEnv * env, jobject paramsJava, Data * d) {
    d->max_iters = getIntUsingGetter(env, paramsJava, "getMaxIters");
    d->eps = getFloatUsingGetter(env, paramsJava, "getEps");
    d->alpha = getFloatUsingGetter(env, paramsJava, "getAlpha");
    d->rho_x = getFloatUsingGetter(env, paramsJava, "getRhoX");
    d->cg_rate = getFloatUsingGetter(env, paramsJava, "getCgRate");
    d->verbose = getBooleanUsingGetter(env, paramsJava, "isVerbose");
    d->normalize = getBooleanUsingGetter(env, paramsJava, "isNormalize");
    d->scale = getFloatUsingGetter(env, paramsJava, "getScale");
    d->warm_start = getBooleanUsingGetter(env, paramsJava, "isWarmStart");
}

Data * getDataStruct(JNIEnv * env, jobject AJava, jdoubleArray bJava, jdoubleArray cJava, jobject paramsJava) {
    Data * d = scs_calloc(1, sizeof(Data));
    d->A = getAMatrix(env, AJava);
    d->b = (*env)->GetDoubleArrayElements(env, bJava, NULL);
    d->m = (*env)->GetArrayLength(env, bJava);
    d->c = (*env)->GetDoubleArrayElements(env, cJava, NULL);
    d->n = (*env)->GetArrayLength(env, cJava);
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

void setSol(JNIEnv * env, jobject solJava, Data * d, Sol * sol) {
    setFloatArrayUsingSetter(env, solJava, sol->x, d->n, "setX");
    setFloatArrayUsingSetter(env, solJava, sol->y, d->m, "setY");
    setFloatArrayUsingSetter(env, solJava, sol->s, d->m, "setS");
}

#ifdef INDIRECTJ
JNIEXPORT void JNICALL Java_scs_IndirectSolver_csolve (JNIEnv *env, jclass clazz, jobject AJava,
        jdoubleArray bJava, jdoubleArray cJava, jobject coneJava, jobject paramsJava, jobject solJava)
#else
JNIEXPORT void JNICALL Java_scs_DirectSolver_csolve (JNIEnv *env, jclass clazz, jobject AJava,
        jdoubleArray bJava, jdoubleArray cJava, jobject coneJava, jobject paramsJava, jobject solJava)
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

    freeData(d,k);
    freeSol(sol);
}

