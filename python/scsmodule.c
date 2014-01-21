#include <Python.h>
#include "scs.h"
#include "numpy/arrayobject.h"

// TODO: when normalizing, make a copy

/* WARNING: this code uses numpy array types
 *
 * WARNING: this code also does not check that the data for the matrix A is
 * actually column compressed storage for a sparse matrix. if it's not, the
 * code will just crash inelegantly. support for cvxopt matrix or scipy sparse
 * is planned, but not likely to be implemented soon.
 */

static PyObject *scsError;
 
static inline int getIntType() {
  switch(NPY_SIZEOF_INT) {
    case 1: return NPY_INT8;
    case 2: return NPY_INT16;
    case 4: return NPY_INT32;
    case 8: return NPY_INT64;
    default: return NPY_INT32;  // defaults to 4 byte int
  }
}

static inline int getDoubleType() {
  return NPY_DOUBLE;
}

static inline void freeDataAndConeOnly(Data *d, Cone *k) {
  // this function is useful since the Data and Cone "structs" do not own the
  // memory for the arrays; numpy does.
  if(d) free(d);
  if(k) free(k);
  d = NULL; k = NULL;
}

// TODO: use PyObject * to keep track of whether or not two objects are equivalent (for warm-starting)
// static const double *prev_Ax, *prev_b, *prev_c;
// static const int *prev_Ai, *prev_Ap, *prev_q;

static Sol *solution = NULL;

static void cleanup()
{
  free_sol(solution);
}

static PyObject *solve(PyObject* self, PyObject *args, PyObject *keywords)
{
  /* Expects a function call 
   *     sol = solve(Ax, Ai, Ap, b, c, f=,l=,q=, MAX_ITERS=, EPS_ABS=, EPS_INFEAS=, ALPHA=, VERBOSE=, NORMALIZE=)
   * The uppercase keywords are optional. If INDIRECT is #define'd, then
   * CG_MAX_ITS and CG_TOL are also optional keyword arguments.
   *
   * "A" is a sparse matrix in column compressed storage. "Ax" are the values,
   * "Ai" are the rows, and "Ap" are the column pointers.
   * `Ax` is a Numpy array of doubles
   * `Ai` is a Numpy array of ints
   * `Ap` is a Numpy array of ints
   *
   * `b` is a (dense) Numpy array of doubles
   * `c` is a (dense) Numpy array of doubles
   *
   * `f` is an integer giving the number of free variables
   * `l` is an integer giving the number of nonnegative constraints
   * `q` is a Numpy array of integers giving the number of cone constraints
   *
   * 
   * MAX_ITERS is an integer. Sets the maximum number of ADMM iterations.
   *  Defaults to 2000.
   * EPS_ABS is a double. Sets the quitting tolerance for ADMM. 
   *  Defaults to 1e-4.
   * EPS_INFEAS is a double. Sets the quitting tolerance for infeasibility.
   *  Defaults to 5e-5.
   * ALPHA is a double in (0,2) (non-inclusive). Sets the over-relaxation
   *  parameter. Defaults to 1.0.
   * VERBOSE is an integer (or Boolean) either 0 or 1. Sets the verbosity of
   *  the solver. Defaults to 1 (or True).
   * NORMALIZE is an integer (or Boolean) either 0 or 1. Tells the solver to
   *  normalize the data. Defaults to 0 (or False).
   *
   * CG_MAX_ITS is an integer. Sets the maximum number of CG iterations.
   *  Defaults to 20.
   * CG_TOL is a double. Sets the tolerance for CG.
   *  Defaults to 1e-3.
   *  
   * The code returns a Python dictionary with three keys, 'x', 'y', and 'status'.
   * These report the primal and dual solution (as numpy arrays) and the solver
   * status (as a string).
   */
     
     
#ifdef INDIRECT
  static char *kwlist[] = {"Ax","Ai","Ap","b","c","f","l","q","MAX_ITERS", "EPS_ABS", "EPS_INFEAS", "ALPHA", "CG_MAX_ITS", "CG_TOL", "VERBOSE", "NORMALIZE", NULL};
#else
  static char *kwlist[] = {"Ax","Ai","Ap","b","c","f","l","q","MAX_ITERS", "EPS_ABS", "EPS_INFEAS", "ALPHA", "VERBOSE", "NORMALIZE", NULL};
#endif
  Data *d = calloc(1,sizeof(Data)); // sets everything to 0
  Cone *k = calloc(1,sizeof(Cone)); // sets everything to 0
  d->MAX_ITERS = 2000;
  d->EPS_ABS = 1e-4;
  d->EPS_INFEAS = 5e-5;
  d->ALPH = 1.0;
  d->VERBOSE = 1;
#ifdef INDIRECT
  d->CG_MAX_ITS = 20;
  d->CG_TOL = 1e-3;
#endif
  
  PyArrayObject *Ax, *Ai, *Ap, *b, *c, *q = NULL;
  PyArrayObject *tmp_arr;
  PyArrayObject *Ax_arr, *Ai_arr, *Ap_arr, *b_arr, *c_arr, *q_arr = NULL;
  
  int intType = getIntType();
  int doubleType = getDoubleType();
  
#ifdef INDIRECT
  if( !PyArg_ParseTupleAndKeywords(args, keywords, "O!O!O!O!O!|iiO!idddidii", kwlist,
      &PyArray_Type, &Ax,
      &PyArray_Type, &Ai,
      &PyArray_Type, &Ap,
      &PyArray_Type, &b,
      &PyArray_Type, &c,
      &(k->f),
      &(k->l),
      &PyArray_Type, &q,
      &(d->MAX_ITERS), 
      &(d->EPS_ABS),
      &(d->EPS_INFEAS),
      &(d->ALPH),
      &(d->CG_MAX_ITS),
      &(d->CG_TOL),
      &(d->VERBOSE),
      &(d->NORMALIZE))
    ) { freeDataAndConeOnly(d,k); return NULL; } 
#else
  if( !PyArg_ParseTupleAndKeywords(args, keywords, "O!O!O!O!O!|iiO!idddii", kwlist,
      &PyArray_Type, &Ax,
      &PyArray_Type, &Ai,
      &PyArray_Type, &Ap,
      &PyArray_Type, &b,
      &PyArray_Type, &c,
      &(k->f),
      &(k->l),
      &PyArray_Type, &q,
      &(d->MAX_ITERS), 
      &(d->EPS_ABS),
      &(d->EPS_INFEAS),
      &(d->ALPH),
      &(d->VERBOSE),
      &(d->NORMALIZE))
    ) { freeDataAndConeOnly(d,k); return NULL; }
#endif
  
  // check that Ax is a double array (and ensure it's of type "double")
  if( !PyArray_ISFLOAT(Ax) ) {
    PyErr_SetString(scsError, "Ax must be a numpy array of floats");
    freeDataAndConeOnly(d,k);
    return NULL;
  }
  tmp_arr = PyArray_GETCONTIGUOUS(Ax);
  Ax_arr = (PyArrayObject *) PyArray_Cast(tmp_arr, doubleType);
  Py_DECREF(tmp_arr);
  d->Ax = (double *) PyArray_DATA(Ax_arr);
  
  // check that Ai is a int array (and ensure it's of type "int")
  if( !PyArray_ISINTEGER(Ai) ) {
    PyErr_SetString(scsError, "Ai must be a numpy array of ints");
    Py_DECREF(Ax_arr);
    freeDataAndConeOnly(d,k);
    return NULL;
  }
  tmp_arr = PyArray_GETCONTIGUOUS(Ai);
  Ai_arr = (PyArrayObject *) PyArray_Cast(tmp_arr, intType);
  Py_DECREF(tmp_arr);
  d->Ai = (int *) PyArray_DATA(Ai_arr);
  
  // check that Ap is a int array (and ensure it's of type "int")
  if( !PyArray_ISINTEGER(Ap) ) {
    PyErr_SetString(scsError, "Ap must be a numpy array of ints");
    Py_DECREF(Ax_arr); Py_DECREF(Ai_arr);
    freeDataAndConeOnly(d,k);
    return NULL;
  }
  tmp_arr = PyArray_GETCONTIGUOUS(Ap);
  Ap_arr = (PyArrayObject *) PyArray_Cast(tmp_arr, intType);
  Py_DECREF(tmp_arr);
  d->Ap = (int *) PyArray_DATA(Ap_arr);
  
  // check that b is a double array (and ensure it's of type "double")
  if( !PyArray_ISFLOAT(b) ) {
    PyErr_SetString(scsError, "b must be a numpy array of floats");
    Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr);
    freeDataAndConeOnly(d,k);
    return NULL;
  }
  tmp_arr = PyArray_GETCONTIGUOUS(b);
  b_arr = (PyArrayObject *) PyArray_Cast(tmp_arr, doubleType);
  Py_DECREF(tmp_arr);
  d->b = (double *) PyArray_DATA(b_arr);
  d->m = (int) PyArray_SIZE(b_arr); // loses precision
  
  // check that c is a double array (and ensure it's of type "double")
  if( !PyArray_ISFLOAT(c) ) {
    PyErr_SetString(scsError, "c must be a numpy array of floats");
    Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr); Py_DECREF(b_arr);
    freeDataAndConeOnly(d,k);
    return NULL;
  }
  tmp_arr = PyArray_GETCONTIGUOUS(c);
  c_arr = (PyArrayObject *) PyArray_Cast(tmp_arr, doubleType);
  Py_DECREF(tmp_arr);
  d->c = (double *) PyArray_DATA(c_arr);
  d->n = (int) PyArray_SIZE(c_arr); // loses precision
  
  // check that q is a int array (and ensure it's of type "double")
  if(q) {
    if( !PyArray_ISINTEGER(q) ) {
      PyErr_SetString(scsError, "q must be a numpy array of ints");
      Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr); Py_DECREF(b_arr); Py_DECREF(c_arr);
      freeDataAndConeOnly(d,k);
      return NULL;
    }
    tmp_arr = PyArray_GETCONTIGUOUS(q);
    q_arr = (PyArrayObject *) PyArray_Cast(tmp_arr, intType);
    Py_DECREF(tmp_arr);
    k->q = (int *) PyArray_DATA(q_arr);
    k->qsize = (int) PyArray_SIZE(q_arr); // loses precision
  }
  
  // TODO: check that parameter values are correct
  
  // solve the problem
  // TODO: preserve the workspace
  solution = scs(d, k);
  
  npy_intp dims[1];
  dims[0] = d->n;
  PyObject *primalSol = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, solution->x);
  dims[0] = d->m;
  PyObject *dualSol = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, solution->y);
  PyObject *returnDict = Py_BuildValue("{s:O,s:O,s:s}","x", primalSol, "y", dualSol, "status", solution->status);
  // give up ownership to the return dictionary
  Py_DECREF(primalSol); Py_DECREF(dualSol); 
  
  // do some cleanup
  freeDataAndConeOnly(d,k);
  Py_DECREF(Ax_arr); Py_DECREF(Ai_arr); Py_DECREF(Ap_arr); Py_DECREF(b_arr); Py_DECREF(c_arr);
  if(q_arr)
    Py_DECREF(q_arr);

  return returnDict;
}

static PyMethodDef scsMethods[] =
{
  {"solve", (PyCFunction)solve, METH_VARARGS | METH_KEYWORDS, 
    "Solve a conic optimization problem."},
  {NULL, NULL, 0, NULL} // sentinel
};

PyMODINIT_FUNC
#ifdef INDIRECT
  initscs_indirect(void)
#else
  initscs_direct(void)
#endif  
{
  PyObject *m;

#ifdef INDIRECT
  m = Py_InitModule("scs_indirect", scsMethods);
#else
  m = Py_InitModule("scs_direct", scsMethods);
#endif
  
  import_array(); // for numpy support
  
  if(m == NULL)
    return;

#ifdef INDIRECT
  scsError = PyErr_NewException("scs_indirect.error", NULL, NULL);
#else
  scsError = PyErr_NewException("scs_direct.error", NULL, NULL);
#endif

  Py_INCREF(scsError);
  PyModule_AddObject(m, "error", scsError);
  
  Py_AtExit(&cleanup);
}
