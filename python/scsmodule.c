#include <Python.h>
#include "glbopts.h"
#include "scs.h"
#include "cones.h"
#include "numpy/arrayobject.h"

/* IMPORTANT: This code now uses numpy array types. It is a private C module
 * in the sense that end users only see the front-facing Python code in
 * "scs.py"; hence, we can get away with the inputs being numpy arrays of
 * the CSR data structures.
 *
 * WARNING: This code also does not check that the data for the sparse 
 * matrices are *actually* in column compressed storage for a sparse matrix. 
 * The C module is not designed to be used stand-alone. If the data provided
 * does not correspond to a CSR matrix, this code will just crash inelegantly.
 * Please use the "solve" interface in scs.py.
 */
//#include "cvxopt.h"

/* Note, Python3.x may require special handling for the idxint and pfloat
 * types. */
static inline int getIntType() {
  switch(sizeof(idxint)) {
    case 1: return NPY_INT8;
    case 2: return NPY_INT16;
    case 4: return NPY_INT32;
    case 8: return NPY_INT64;
    default: return NPY_INT32;  // defaults to 4 byte int
  }
}

static inline int getDoubleType() {
  // known bug, if pfloat isn't "double", will cause aliasing in memory
  return NPY_DOUBLE;
}

static inline PyArrayObject *getContiguous(PyArrayObject *array, int typenum) {
  // gets the pointer to the block of contiguous C memory
  // the overhead should be small unless the numpy array has been
  // reordered in some way or the data type doesn't quite match
  //
  // the "new_owner" pointer has to have Py_DECREF called on it; it owns
  // the "new" array object created by PyArray_Cast
  //
  static PyArrayObject *tmp_arr;
  PyArrayObject *new_owner;
  tmp_arr = PyArray_GETCONTIGUOUS(array);
  new_owner = (PyArrayObject *) PyArray_Cast(tmp_arr, typenum);
  Py_DECREF(tmp_arr);
  return new_owner;
}

static int printErr(char * key) {
    char str[80];
    sprintf(str, "error parsing '%s'", key);
    PyErr_SetString(PyExc_TypeError, str);
    return -1; 
}

static int getConeIntDim(char * key, idxint * v, PyObject * cone) {
    /* get cone['l'] */
    *v = 0;
    PyObject *obj = PyDict_GetItemString(cone, key);
    if(obj) {
        if(!PyInt_Check(obj) || !((*v = (idxint) PyInt_AsLong(obj)) >= 0)) {
            return printErr(key);
        }
    }
    return 0;
}

static int getConeArrDim(char * key, idxint ** varr, idxint * vsize, PyObject * cone) {
    /* get cone['key'] */
    idxint i, n = 0;
    idxint * q = NULL;
    PyObject *obj = PyDict_GetItemString(cone, key);
    if(obj) {
        if(PyList_Check(obj)) {
            n = PyList_Size(obj);
            q = calloc(n, sizeof(idxint));
            for (i = 0; i < n; ++i) {
                PyObject *qi = PyList_GetItem(obj, i);
                if (!PyInt_Check(qi) || !((q[i] = (idxint) PyInt_AsLong(qi)) > 0)) {
                    return printErr(key);
                }
            }
        } else if (PyInt_Check(obj)) {
            n = 1;
            q = malloc(sizeof(idxint));
            if (!PyInt_Check(obj) || !((*q = (idxint) PyInt_AsLong(obj)) > 0)) {
                return printErr(key);
            }
        } else {
            return printErr(key);
        }
    } 
    *vsize = n;
    *varr = q;
    return 0;
}

static int getOptIntParam(char * key, idxint * v, idxint defVal, PyObject * opts){
    *v = defVal;
    if (opts) {
        PyObject *obj = PyDict_GetItemString(opts, key);
        if (obj) {
            if(!PyInt_Check(obj) || !((*v = (idxint) PyInt_AsLong(obj)) >= 0)) {
                return printErr(key);
            }
        }
    }
    return 0;
}

static int getOptFloatParam(char * key, pfloat * v, pfloat defVal, PyObject * opts){
    *v = defVal;
    if (opts) {
        PyObject *obj = PyDict_GetItemString(opts, key);
        if (obj) {
            if(!PyFloat_Check(obj) || !((*v = (pfloat) PyFloat_AsDouble(obj)) >= 0)) {
                char str[80];
                sprintf(str, "'%s' ought to be a nonnegative float", key);
                PyErr_SetString(PyExc_TypeError, str);
                return -1;
            }
        }
    }
    return 0;
}

static int parseOpts(Data *d, PyObject * opts) {
    if (getOptIntParam("MAX_ITERS", &(d->MAX_ITERS), 2000, opts) < 0)
        return -1;
    if (getOptIntParam("VERBOSE", &(d->VERBOSE), 1, opts) < 0)
        return -1;
    if (getOptIntParam("NORMALIZE", &(d->NORMALIZE), 1, opts))
        return -1;
    if (getOptFloatParam("EPS", &(d->EPS), 1e-3, opts) < 0)
        return -1;
    if (getOptFloatParam("ALPHA", &(d->ALPHA), 1.8, opts) < 0)
        return -1;
    if (getOptFloatParam("UNDET_TOL", &(d->UNDET_TOL), 1e-9, opts) < 0)
        return -1;
    if (getOptFloatParam("RHO_X", &(d->RHO_X), 1e-3, opts) < 0)
        return -1;
    return 0;
}

/* The PyInt variable is a PyLong in Python3.x.
 */
#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong PyLong_AsLong
#define PyInt_Check PyLong_Check
#endif

static PyObject *csolve(PyObject* self, PyObject *args, PyObject *kwargs)
{
  /* Expects a function call
   *     sol = csolve((m,n,p),c,Gx,Gi,Gp,h,dims,Ax,Ai,Ap,b,verbose)
   * where
   *
   * the triple (m,n,p) corresponds to:
   *    `m`: the rows of G
   *    `n`: the cols of G and A, must agree with the length of c
   *    `p`: the rows of A
   * `c` is a Numpy array of pfloats
   * "G" is a sparse matrix in column compressed storage. "Gx" are the values,
   * "Gi" are the rows, and "Gp" are the column pointers.
   * `Gx` is a Numpy array of pfloats
   * `Gi` is a Numpy array of ints
   * `Gp` is a Numpy array of ints
   * `h` is a Numpy array
   * `dims` is a dictionary with
   *    `dims['l']` an integer specifying the dimension of positive orthant cone
   *    `dims['q']` an *list* specifying dimensions of second-order cones
   *
   * "A" is an optional sparse matrix in column compressed storage. "Ax" are 
   * the values, "Ai" are the rows, and "Ap" are the column pointers.
   * `Ax` is a Numpy array of pfloats
   * `Ai` is a Numpy array of ints
   * `Ap` is a Numpy array of ints
   * `b` is an optional argument, which is a Numpy array of pfloats
   * `verbose` is an optional bool signaling whether to print info
   *
   * This call will solve the problem
   *
   *    minimize     c'*x
   *    subject to   A*x = b
   *                 h - G*x \in K
   *
   * The code returns a Python dictionary with five keys, 'x', 'y', 'info', 's',
   * and 'z'. These correspond to the following:
   *
   * `x`: primal variables
   * `y`: dual variables for equality constraints
   * `s`: slacks for Gx + s <= h, s \in K
   * `z`: dual variables for inequality constraints s \in K
   * `info`: another dictionary with the following fields:
   *    exitflag: 0=OPTIMAL, 1=PRIMAL INFEASIBLE, 2=DUAL INFEASIBLE, -1=MAXIT REACHED
   *  infostring: gives information about the status of solution
   *       pcost: value of primal objective
   *       dcost: value of dual objective
   *        pres: primal residual on inequalities and equalities
   *        dres: dual residual
   *        pinf: primal infeasibility measure
   *        dinf: dual infeasibility measure
   *     pinfres: NaN
   *     dinfres: 3.9666e+15
   *         gap: duality gap
   *      relgap: relative duality gap
   *          r0: ???
   *      numerr: numerical error?
   *        iter: number of iterations
   *      timing: dictionary with timing information
   */

  /* data structures for arguments */
  //matrix *c, *h, *b = NULL;
  //spmatrix *G, *A = NULL;
  
  PyArrayObject *Ax, *Ai, *Ap, *c, *b;
  PyObject *cone, *opts;

  /* scs data structures */
  
  Data * d = calloc(sizeof(Data),1); 
  Cone * k = calloc(sizeof(Cone),1); 

  static char *kwlist[] = {"shape", "Ax", "Ai", "Ap", "b", "c", "cone", "opts", NULL};
  // parse the arguments and ensure they are the correct type
#ifdef DLONG
  static char *argparse_string = "(ll)O!O!O!O!O!O!|O!";
#else
  static char *argparse_string = "(ii)O!O!O!O!O!O!|O!";
#endif
    
  if( !PyArg_ParseTupleAndKeywords(args, kwargs, argparse_string, kwlist,
      &(d->m), &(d->n),
      &PyArray_Type, &Ax,
      &PyArray_Type, &Ai,
      &PyArray_Type, &Ap,
      &PyArray_Type, &b,
      &PyArray_Type, &c,
      &PyDict_Type, &cone,
      &PyDict_Type, &opts)
    ) { return NULL; }
  
  if (d->m < 0) {
    PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
    return NULL;
  }

  if (d->n < 0) {
    PyErr_SetString(PyExc_ValueError, "n must be a positive integer");
    return NULL;
  }
  
  /* get the typenum for the primitive idxint and pfloat types */
  int intType = getIntType();
  int pfloatType = getDoubleType();

  /* set G */
  if( !PyArray_ISFLOAT(Ax) || PyArray_NDIM(Ax) != 1) {
    PyErr_SetString(PyExc_TypeError, "Ax must be a numpy array of floats");
    return NULL;
  }
  if( !PyArray_ISINTEGER(Ai) || PyArray_NDIM(Ai) != 1) {
    PyErr_SetString(PyExc_TypeError, "Ai must be a numpy array of ints");
    return NULL;
  }
  if( !PyArray_ISINTEGER(Ap) || PyArray_NDIM(Ap) != 1) {
    PyErr_SetString(PyExc_TypeError, "Ap must be a numpy array of ints");
    return NULL;
  }
  PyArrayObject *Ax_arr = getContiguous(Ax, pfloatType);
  PyArrayObject *Ai_arr = getContiguous(Ai, intType);
  PyArrayObject *Ap_arr = getContiguous(Ap, intType);
  d->Ax = (pfloat *) PyArray_DATA(Ax_arr);
  d->Ai = (idxint *) PyArray_DATA(Ai_arr);
  d->Ap = (idxint *) PyArray_DATA(Ap_arr);
  //d->Anz = d->Ap[d->n];
  d->Anz = PyArray_DIM(Ai,0);
  /* set c */
  if (!PyArray_ISFLOAT(c) || PyArray_NDIM(c) != 1) {
      PyErr_SetString(PyExc_TypeError, "c must be a dense numpy array with one dimension");
      Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
      return NULL;
  }
  if (PyArray_DIM(c,0) != d->n){
      PyErr_SetString(PyExc_ValueError, "c has incompatible dimension with A");
      Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
      return NULL;
  }
  PyArrayObject *c_arr = getContiguous(c, pfloatType);
  d->c = (pfloat *) PyArray_DATA(c_arr);

  /* set b */
  if (!PyArray_ISFLOAT(b) || PyArray_NDIM(b) != 1) {
      PyErr_SetString(PyExc_TypeError, "b must be a dense numpy array with one dimension");
      Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
      Py_DECREF(c_arr);
      return NULL;
  }
  if (PyArray_DIM(b,0) != d->m){
      PyErr_SetString(PyExc_ValueError, "b has incompatible dimension with A");
      Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
      Py_DECREF(c_arr);
      return NULL;
  }
  PyArrayObject *b_arr = getContiguous(b, pfloatType);
  d->b = (pfloat *) PyArray_DATA(b_arr);

  getConeIntDim("f", &(k->f), cone);
  getConeIntDim("l", &(k->l), cone);
  getConeArrDim("q", &(k->q), &(k->qsize), cone);
  getConeArrDim("s", &(k->s), &(k->ssize), cone);
  getConeIntDim("ep", &(k->ep), cone);
  getConeIntDim("ed", &(k->ed), cone);
  parseOpts(d,opts);
  
  /* Solve! */
  Sol sol;
  Info info;
  scs(d, k, &sol, &info);

  /* create output (all data is *deep copied*) */
  // TODO: request CVXOPT API for constructing from existing pointer
  /* x */
  // matrix *x;
  // if(!(x = Matrix_New(n,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(x), mywork->x, n*sizeof(pfloat));
  npy_intp veclen[1];
  veclen[0] = d->n;
  PyObject *x = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.x);

  /* y */
  // matrix *y;
  // if(!(y = Matrix_New(p,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(y), mywork->y, p*sizeof(pfloat));
  veclen[0] = d->m;
  PyObject *y = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.y);

  /* s */
  // matrix *s;
  // if(!(s = Matrix_New(m,1,DOUBLE)))
  //   return PyErr_NoMemory();
  // memcpy(MAT_BUFD(s), mywork->s, m*sizeof(pfloat));
  veclen[0] = d->m;
  PyObject *s = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.s);
  
    PyObject *infoDict = Py_BuildValue(
    "{s:l,s:l,s:d,s:d,s:d,s:d,s:d,s:d,s:s}",
    "statusVal", (idxint)info.statusVal,
    "iter", (idxint)info.iter,
    "pobj", (pfloat)info.pobj,
    "dobj", (pfloat)info.dobj,
    "resPri", (pfloat)info.resPri,
    "resDual", (pfloat)info.resDual,
    "relGap", (pfloat)info.relGap,
    "time",(pfloat)(info.time/1e3),
    "status",info.status);

  PyObject *returnDict = Py_BuildValue(
    "{s:O,s:O,s:O,s:O}",
    "x",x,
    "y",y,
    "s",s,
    "info",infoDict);
  // give up ownership to the return dictionary
  Py_DECREF(x); Py_DECREF(y); Py_DECREF(s); Py_DECREF(infoDict);

  // no longer need pointers to arrays that held primitives
  if(k) {
      if(k->q) scs_free(k->q);
      if(k->s) scs_free(k->s);
      scs_free(k);
  }   
  Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
  Py_DECREF(c_arr); Py_DECREF(b_arr);
  scs_free(d);
  return returnDict;
}

static PyMethodDef scsMethods[] =
{
  {"csolve", (PyCFunction)csolve, METH_VARARGS | METH_KEYWORDS,
    "Solve a convex cone problem using scs."},
  {NULL, NULL, 0, NULL} // sentinel
};

/* Module initialization */
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef = {
          PyModuleDef_HEAD_INIT,
          "_scs",             /* m_name */
          "Solve a convex cone problem using scs.",  /* m_doc */
          -1,                  /* m_size */
          scsMethods,         /* m_methods */
          NULL,                /* m_reload */
          NULL,                /* m_traverse */
          NULL,                /* m_clear */
          NULL,                /* m_free */
  };
#endif

static PyObject* moduleinit(void)
{
  PyObject* m;

#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
#else
  #ifdef INDIRECT
    m = Py_InitModule("_scs_indirect", scsMethods);
  #else
    m = Py_InitModule("_scs_direct", scsMethods);
  #endif
#endif
  
  //if (import_array() < 0) return NULL; // for numpy arrays
  //if (import_cvxopt() < 0) return NULL; // for cvxopt support

  if (m == NULL)
    return NULL;

  return m;
};

#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC
  #ifdef INDIRECT
    PyInit__scs_indirect(void)
  #else
    PyInit__scs_direct(void)
  #endif
  {
    import_array(); // for numpy arrays
    return moduleinit();
  }
#else
  PyMODINIT_FUNC 
  #ifdef INDIRECT
    init_scs_indirect(void)
  #else
    init_scs_direct(void)
  #endif
  {
    import_array(); // for numpy arrays
    moduleinit();
  }
#endif
