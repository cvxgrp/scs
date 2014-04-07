#include <Python.h>
#include "glbopts.h"
#include "scs.h"
#include "cones.h"
#include "../linsys/amatrix.h"
#include "numpy/arrayobject.h"

/* IMPORTANT: This code now uses numpy array types. It is a private C module
 * in the sense that end users only see the front-facing Python code in
 * "scs.py"; hence, we can get away with the inputs being numpy arrays of
 * the CSC data structures.
 *
 * WARNING: This code also does not check that the data for the sparse
 * matrices are *actually* in column compressed storage for a sparse matrix.
 * The C module is not designed to be used stand-alone. If the data provided
 * does not correspond to a CSR matrix, this code will just crash inelegantly.
 * Please use the "solve" interface in scs.py.
 */

static int intType;
static int pfloatType;

struct ScsPyData {
	PyArrayObject * Ax;
	PyArrayObject * Ai;
	PyArrayObject * Ap;
	PyArrayObject * b;
	PyArrayObject * c;
	PyArrayObject * x0;
	PyArrayObject * y0;
	PyArrayObject * s0;
};

/* Note, Python3.x may require special handling for the idxint and pfloat
 * types. */
static int getIntType(void) {
	switch (sizeof(idxint)) {
	case 1:
		return NPY_INT8;
	case 2:
		return NPY_INT16;
	case 4:
		return NPY_INT32;
	case 8:
		return NPY_INT64;
	default:
		return NPY_INT32; /* defaults to 4 byte int */
	}
}

static int getDoubleType(void) {
	/* known bug, if pfloat isn't "double", will cause aliasing in memory */
	return NPY_DOUBLE;
}

static PyArrayObject *getContiguous(PyArrayObject *array, int typenum) {
	/* gets the pointer to the block of contiguous C memory */
	/* the overhead should be small unless the numpy array has been */
	/* reordered in some way or the data type doesn't quite match */
	/* */
	/* the "new_owner" pointer has to have Py_DECREF called on it; it owns */
	/* the "new" array object created by PyArray_Cast */
	/* */
	static PyArrayObject *tmp_arr;
	PyArrayObject *new_owner;
	tmp_arr = PyArray_GETCONTIGUOUS(array);
	new_owner = (PyArrayObject *) PyArray_Cast(tmp_arr, typenum);
	Py_DECREF(tmp_arr);
	return new_owner;
}

static int printErr(char * key) {
	scs_printf("error parsing '%s'", key);
	return -1;
}

static idxint getWarmStart(char * key, pfloat ** x, PyArrayObject ** px0, idxint l, PyObject * warm) {
	*x = scs_calloc(l, sizeof(pfloat));
	PyArrayObject *x0 = (PyArrayObject *) PyDict_GetItemString(warm, key);
	if (x0) {
		if (!PyArray_ISFLOAT(x0) || PyArray_NDIM(x0) != 1 || PyArray_DIM(x0,0) != l) {
			scs_printf("Error parsing warm-start input\n");
			return 0;
		} else {
			*px0 = getContiguous(x0, pfloatType);
			*x = (pfloat *) PyArray_DATA(*px0);
			return 1;
		}
	}
	return 0;
}

static int getConeIntDim(char * key, idxint * v, PyObject * cone) {
	/* get cone['l'] */
	*v = 0;
	PyObject *obj = PyDict_GetItemString(cone, key);
	if (obj) {
		if (!PyInt_Check(obj) || !((*v = (idxint) PyInt_AsLong(obj)) >= 0)) {
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
	if (obj) {
		if (PyList_Check(obj)) {
			n = PyList_Size(obj);
			q = scs_calloc(n, sizeof(idxint));
			for (i = 0; i < n; ++i) {
				PyObject *qi = PyList_GetItem(obj, i);
				if (!PyInt_Check(qi) || !((q[i] = (idxint) PyInt_AsLong(qi)) >= 0)) {
					return printErr(key);
				}
			}
		} else if (PyInt_Check(obj)) {
			n = 1;
			q = scs_malloc(sizeof(idxint));
			if (!PyInt_Check(obj) || !((*q = (idxint) PyInt_AsLong(obj)) >= 0)) {
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

static int getOptIntParam(char * key, idxint * v, idxint defVal, PyObject * opts) {
	*v = defVal;
	if (opts) {
		PyObject *obj = PyDict_GetItemString(opts, key);
		if (obj) {
			if (!PyInt_Check(obj) || !((*v = (idxint) PyInt_AsLong(obj)) >= 0)) {
				return printErr(key);
			}
		}
	}
	return 0;
}

static int getOptFloatParam(char * key, pfloat * v, pfloat defVal, PyObject * opts) {
	*v = defVal;
	if (opts) {
		PyObject *obj = PyDict_GetItemString(opts, key);
		if (obj) {
			if (PyInt_Check(obj)) {
				if ((*v = (pfloat) PyInt_AsLong(obj)) < 0) {
					scs_printf("ERROR: '%s' ought to be a nonnegative float\n", key);
					return -1;
				}
			} else if (!PyFloat_Check(obj) || !((*v = (pfloat) PyFloat_AsDouble(obj)) >= 0)) {
				scs_printf("ERROR: '%s' ought to be a nonnegative float\n", key);
				return -1;
			}
		}
	}
	return 0;
}

static int parseOpts(Data *d, PyObject * opts) {
	if (getOptIntParam("MAX_ITERS", &(d->MAX_ITERS), 2500, opts) < 0)
		return -1;
	if (getOptIntParam("VERBOSE", &(d->VERBOSE), 1, opts) < 0)
		return -1;
	if (getOptIntParam("NORMALIZE", &(d->NORMALIZE), 1, opts))
		return -1;
	if (getOptFloatParam("SCALE", &(d->SCALE), 5, opts) < 0)
		return -1;
	if (getOptFloatParam("EPS", &(d->EPS), 1e-3, opts) < 0)
		return -1;
	if (getOptFloatParam("CG_RATE", &(d->CG_RATE), 1.5, opts) < 0)
		return -1;
	if (getOptFloatParam("ALPHA", &(d->ALPHA), 1.8, opts) < 0)
		return -1;
	if (getOptFloatParam("RHO_X", &(d->RHO_X), 1e-3, opts) < 0)
		return -1;
	return 0;
}

static void freePyData(Data * d, Cone * k, struct ScsPyData * ps) {
	if (ps->Ax)
		Py_DECREF(ps->Ax);
	if (ps->Ai)
		Py_DECREF(ps->Ai);
	if (ps->Ap)
		Py_DECREF(ps->Ap);
	if (ps->b)
		Py_DECREF(ps->b);
	if (ps->c)
		Py_DECREF(ps->c);
	if (ps->x0)
		Py_DECREF(ps->x0);
	if (ps->y0)
		Py_DECREF(ps->y0);
	if (ps->s0)
		Py_DECREF(ps->s0);
	if (k) {
		if (k->q)
			scs_free(k->q);
		if (k->s)
			scs_free(k->s);
		scs_free(k);
	}
	if (d) {
		if (d->A)
			scs_free(d->A);
		scs_free(d);
	}
}

static PyObject * finishWithErr(Data * d, Cone * k, struct ScsPyData * ps, char * str) {
	PyErr_SetString(PyExc_ValueError, str);
	freePyData(d, k, ps);
	return NULL;
}

/* The PyInt variable is a PyLong in Python3.x.
 */
#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong PyLong_AsLong
#define PyInt_Check PyLong_Check
#endif
static PyObject *csolve(PyObject* self, PyObject *args, PyObject *kwargs) {
	/* Expects a function call
	 *     sol = csolve((m,n),Ax,Ai,Ap,b,c,cone,opts)
	 * where
	 *
	 * the pair (m,n) corresponds to:
	 *    `m`: the rows of A, must agree with the length of b
	 *    `n`: the cols of A, must agree with the length of c
	 * `c` is a Numpy array of doubles
	 * "A" is a sparse matrix in column compressed storage. "Ax" are the values,
	 * "Ai" are the rows, and "Ap" are the column pointers.
	 * `Ax` is a Numpy array of doubles
	 * `Ai` is a Numpy array of ints
	 * `Ap` is a Numpy array of ints
	 * `b` is a Numpy array of doubles
	 * `cone` is a dictionary with
	 *    `cone['l']` an integer specifying the dimension of positive orthant cone
	 *    `cone['q']` a *list* specifying dimensions of second-order cones
	 *
	 * This call will solve the problem
	 *
	 *    minimize     c'*x
	 *    subject to   A*x + s = b
	 *                 s \in K
	 *
	 * The code returns a Python dictionary with keys 'x', 'y', 's' and 'info'.
	 * These correspond to the following:
	 *
	 * `x`: primal variables
	 * `y`: dual variables for conic constraints
	 * `s`: slacks for Ax + s = b, s \in K
	 * `info`: another dictionary with the following fields:
	 */

	/* data structures for arguments */
	PyArrayObject *Ax, *Ai, *Ap, *c, *b;
	PyObject *cone, *opts, *warm = NULL;
	struct ScsPyData ps = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
	/* scs data structures */
	Data * d = scs_calloc(sizeof(Data), 1);
	Cone * k = scs_calloc(sizeof(Cone), 1);

	static char *kwlist[] = { "shape", "Ax", "Ai", "Ap", "b", "c", "cone", "opts", "warm", NULL };
	/* parse the arguments and ensure they are the correct type */
#ifdef DLONG
	static char *argparse_string = "(ll)O!O!O!O!O!O!|O!O!";
#else
	static char *argparse_string = "(ii)O!O!O!O!O!O!|O!O!";
#endif

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, argparse_string, kwlist, &(d->m), &(d->n), &PyArray_Type, &Ax,
			&PyArray_Type, &Ai, &PyArray_Type, &Ap, &PyArray_Type, &b, &PyArray_Type, &c, &PyDict_Type, &cone,
			&PyDict_Type, &opts, &PyDict_Type, &warm)) {
		return NULL;
	}

	if (d->m < 0) {
		PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
		return NULL;
	}

	if (d->n < 0) {
		PyErr_SetString(PyExc_ValueError, "n must be a positive integer");
		return NULL;
	}

	/* get the typenum for the primitive idxint and pfloat types */
	intType = getIntType();
	pfloatType = getDoubleType();

	/* set A */
	if (!PyArray_ISFLOAT(Ax) || PyArray_NDIM(Ax) != 1) {
		return finishWithErr(d, k, &ps, "Ax must be a numpy array of floats");
	}
	if (!PyArray_ISINTEGER(Ai) || PyArray_NDIM(Ai) != 1) {
		return finishWithErr(d, k, &ps, "Ai must be a numpy array of ints");
	}
	if (!PyArray_ISINTEGER(Ap) || PyArray_NDIM(Ap) != 1) {
		return finishWithErr(d, k, &ps, "Ap must be a numpy array of ints");
	}
	ps.Ax = getContiguous(Ax, pfloatType);
	ps.Ai = getContiguous(Ai, intType);
	ps.Ap = getContiguous(Ap, intType);
	AMatrix * A = scs_malloc(sizeof(AMatrix));
	A->x = (pfloat *) PyArray_DATA(ps.Ax);
	A->i = (idxint *) PyArray_DATA(ps.Ai);
	A->p = (idxint *) PyArray_DATA(ps.Ap);
	d->A = A;
	/*d->Anz = d->Ap[d->n]; */
	/*d->Anz = PyArray_DIM(Ai,0); */
	/* set c */
	if (!PyArray_ISFLOAT(c) || PyArray_NDIM(c) != 1) {
		return finishWithErr(d, k, &ps, "c must be a dense numpy array with one dimension");
	}
	if (PyArray_DIM(c,0) != d->n) {
		return finishWithErr(d, k, &ps, "c has incompatible dimension with A");
	}
	ps.c = getContiguous(c, pfloatType);
	d->c = (pfloat *) PyArray_DATA(ps.c);
	/* set b */
	if (!PyArray_ISFLOAT(b) || PyArray_NDIM(b) != 1) {
		return finishWithErr(d, k, &ps, "b must be a dense numpy array with one dimension");
	}
	if (PyArray_DIM(b,0) != d->m) {
		return finishWithErr(d, k, &ps, "b has incompatible dimension with A");
	}
	ps.b = getContiguous(b, pfloatType);
	d->b = (pfloat *) PyArray_DATA(ps.b);

	if (getConeIntDim("f", &(k->f), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field f");
	}
	if (getConeIntDim("l", &(k->l), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field l");
	}
	if (getConeArrDim("q", &(k->q), &(k->qsize), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field q");
	}
	if (getConeArrDim("s", &(k->s), &(k->ssize), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field s");
	}
	if (getConeIntDim("ep", &(k->ep), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field ep");
	}
	if (getConeIntDim("ed", &(k->ed), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field ed");
	}
	if (parseOpts(d, opts) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse opts");
	}

	/* Solve! */
	Sol sol = { 0 };
	d->WARM_START = 0;
	if (warm) {
		d->WARM_START = getWarmStart("x", &(sol.x), &(ps.x0), d->n, warm);
		d->WARM_START |= getWarmStart("y", &(sol.y), &(ps.y0), d->m, warm);
		d->WARM_START |= getWarmStart("s", &(sol.s), &(ps.s0), d->m, warm);
	}
	Info info;
	scs(d, k, &sol, &info);

	/* create output (all data is *deep copied*) */
	/* x */
	/* matrix *x; */
	/* if(!(x = Matrix_New(n,1,DOUBLE))) */
	/*   return PyErr_NoMemory(); */
	/* memcpy(MAT_BUFD(x), mywork->x, n*sizeof(pfloat)); */
	npy_intp veclen[1];
	veclen[0] = d->n;
	PyObject *x = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.x);

	/* y */
	/* matrix *y; */
	/* if(!(y = Matrix_New(p,1,DOUBLE))) */
	/*   return PyErr_NoMemory(); */
	/* memcpy(MAT_BUFD(y), mywork->y, p*sizeof(pfloat)); */
	veclen[0] = d->m;
	PyObject *y = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.y);

	/* s */
	/* matrix *s; */
	/* if(!(s = Matrix_New(m,1,DOUBLE))) */
	/*   return PyErr_NoMemory(); */
	/* memcpy(MAT_BUFD(s), mywork->s, m*sizeof(pfloat)); */
	veclen[0] = d->m;
	PyObject *s = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.s);

	PyObject *infoDict = Py_BuildValue("{s:l,s:l,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:s}", "statusVal",
			(idxint) info.statusVal, "iter", (idxint) info.iter, "pobj", (pfloat) info.pobj, "dobj", (pfloat) info.dobj,
			"resPri", (pfloat) info.resPri, "resDual", (pfloat) info.resDual, "relGap", (pfloat) info.relGap,
			"solveTime", (pfloat) (info.solveTime / 1e3), "setupTime", (pfloat) (info.setupTime / 1e3), "status",
			info.status);

	PyObject *returnDict = Py_BuildValue("{s:O,s:O,s:O,s:O}", "x", x, "y", y, "s", s, "info", infoDict);
	/* give up ownership to the return dictionary */
	Py_DECREF(x);
	Py_DECREF(y);
	Py_DECREF(s);
	Py_DECREF(infoDict);

	/* no longer need pointers to arrays that held primitives */
	freePyData(d, k, &ps);
	return returnDict;
}

static PyMethodDef scsMethods[] = { { "csolve", (PyCFunction) csolve, METH_VARARGS | METH_KEYWORDS,
		"Solve a convex cone problem using scs." }, { NULL, NULL, 0, NULL } /* sentinel */
};

/* Module initialization */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"_scs", /* m_name */
	"Solve a convex cone problem using scs.", /* m_doc */
	-1, /* m_size */
	scsMethods, /* m_methods */
	NULL, /* m_reload */
	NULL, /* m_traverse */
	NULL, /* m_clear */
	NULL, /* m_free */
};
#endif

static PyObject* moduleinit(void) {
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

	/*if (import_array() < 0) return NULL; // for numpy arrays */
	/*if (import_cvxopt() < 0) return NULL; // for cvxopt support */

	if (m == NULL)
		return NULL;

	return m;
}
;

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC
#ifdef INDIRECT
PyInit__scs_indirect(void)
#else
PyInit__scs_direct(void)
#endif
{
	import_array(); /* for numpy arrays */
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
	import_array(); /* for numpy arrays */
	moduleinit();
}
#endif
