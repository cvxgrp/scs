#include <Python.h>
#include "glbopts.h"
#include "scs.h"
#include "cones.h"
#include "linsys/amatrix.h"
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

/* The PyInt variable is a PyLong in Python3.x.
 */
#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong PyLong_AsLong
#define PyInt_Check PyLong_Check
#endif


static int intType;
static int scs_floatType;

struct ScsPyData {
	PyArrayObject * Ax;
	PyArrayObject * Ai;
	PyArrayObject * Ap;
	PyArrayObject * b;
	PyArrayObject * c;
};

/* Note, Python3.x may require special handling for the scs_int and scs_float
 * types. */
static int getIntType(void) {
	switch (sizeof(scs_int)) {
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
	/* known bug, if scs_float isn't "double", will cause aliasing in memory */
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
	PySys_WriteStderr("error parsing '%s'\n", key);
	return -1;
}

/* gets warm starts from warm dict, doesn't destroy input warm start data */
static scs_int getWarmStart(char * key, scs_float ** x, scs_int l, PyObject * warm) {
    PyArrayObject *x0 = (PyArrayObject *) PyDict_GetItemString(warm, key);
	*x = scs_calloc(l, sizeof(scs_float));
    if (x0) {
		if (!PyArray_ISFLOAT(x0) || PyArray_NDIM(x0) != 1 || PyArray_DIM(x0, 0) != l) {
			PySys_WriteStderr("Error parsing warm-start input\n");
			return 0;
		} else {
			PyArrayObject * px0 = getContiguous(x0, scs_floatType);
			memcpy(*x, (scs_float *) PyArray_DATA(px0), l * sizeof(scs_float));
            Py_DECREF(px0);
			return 1;
		}
	}
	return 0;
}

static int getConeArrDim(char * key, scs_int ** varr, scs_int * vsize, PyObject * cone) {
	/* get cone['key'] */
	scs_int i, n = 0;
	scs_int * q = NULL;
	PyObject *obj = PyDict_GetItemString(cone, key);
	if (obj) {
		if (PyList_Check(obj)) {
			n = PyList_Size(obj);
			q = scs_calloc(n, sizeof(scs_int));
			for (i = 0; i < n; ++i) {
				PyObject *qi = PyList_GetItem(obj, i);
	            if ( !( (PyInt_Check(qi) && ((q[i] = (scs_int) PyInt_AsLong(qi)) >= 0)) ||
	                        (PyLong_Check(qi) && ((q[i] = (scs_int) PyLong_AsLong(qi)) >= 0))) ) {
					return printErr(key);
	            }
			}
		} else if (PyInt_Check(obj) || PyLong_Check(obj)) {
			n = 1;
			q = scs_malloc(sizeof(scs_int));
            if ( !( (PyInt_Check(obj) && ((*q = (scs_int) PyInt_AsLong(obj)) >= 0)) ||
                        (PyLong_Check(obj) && ((*q = (scs_int) PyLong_AsLong(obj)) >= 0))) ) {
				return printErr(key);
			}
		} else {
			return printErr(key);
        }
        if (PyErr_Occurred()) {
            /* potentially could have been triggered before */
            return printErr(key);
        }
    }
    *vsize = n;
    *varr = q;
    return 0;
}

static int getConeFloatArr(char * key, scs_float ** varr, scs_int * vsize, PyObject * cone) {
    /* get cone['key'] */
    scs_int i, n = 0;
    scs_float * q = NULL;
    PyObject *obj = PyDict_GetItemString(cone, key);
    if (obj) {
        if (PyList_Check(obj)) {
            n = PyList_Size(obj);
            q = scs_calloc(n, sizeof(scs_float));
            for (i = 0; i < n; ++i) {
                PyObject *qi = PyList_GetItem(obj, i);
                q[i] = (scs_float) PyFloat_AsDouble(qi);
            }
        } else if (PyInt_Check(obj) || PyLong_Check(obj) || PyFloat_Check(obj)) {
            n = 1;
            q = scs_malloc(sizeof(scs_float));
            q[0] = (scs_float) PyFloat_AsDouble(obj);
        } else {
            return printErr(key);
        }
        if (PyErr_Occurred()) {
            /* potentially could have been triggered before */
            return printErr(key);
        }
    }
    *vsize = n;
    *varr = q;
    return 0;
}


static int getPosIntParam(char * key, scs_int * v, scs_int defVal, PyObject * opts) {
    *v = defVal;
    if (opts) {
        PyObject *obj = PyDict_GetItemString(opts, key);
        if (obj) {
            if ( !( (PyInt_Check(obj) && ((*v = (scs_int) PyInt_AsLong(obj)) >= 0)) ||
                        (PyLong_Check(obj) && ((*v = (scs_int) PyLong_AsLong(obj)) >= 0))) ) {
                return printErr(key);
            }
        }
    }
    return 0;
}

static void freePyData(Data * d, Cone * k, struct ScsPyData * ps) {
    if (ps->Ax) {
        Py_DECREF(ps->Ax);
    }
    if (ps->Ai) {
        Py_DECREF(ps->Ai);
    }
    if (ps->Ap) {
        Py_DECREF(ps->Ap);
    }
    if (ps->b) {
        Py_DECREF(ps->b);
    }
    if (ps->c) {
        Py_DECREF(ps->c);
    }
    if (k) {
        if (k->q)
            scs_free(k->q);
        if (k->s)
            scs_free(k->s);
        if (k->p)
            scs_free(k->p);
        scs_free(k);
    }
    if (d) {
        if (d->A)
            scs_free(d->A);
        if (d->stgs)
            scs_free(d->stgs);
        scs_free(d);
    }
}

static PyObject * finishWithErr(Data * d, Cone * k, struct ScsPyData * ps, char * str) {
	PyErr_SetString(PyExc_ValueError, str);
	freePyData(d, k, ps);
	return NULL;
}

static PyObject *version(PyObject* self) {
  return Py_BuildValue("s", scs_version());
}

static PyObject *csolve(PyObject* self, PyObject *args, PyObject *kwargs) {
	/* data structures for arguments */
	PyArrayObject *Ax, *Ai, *Ap, *c, *b;
	PyObject *cone, *warm = NULL;
    PyObject *verbose = NULL;
    PyObject *normalize = NULL;
	struct ScsPyData ps = { NULL, NULL, NULL, NULL, NULL, };
	/* scs data structures */
	Data * d = scs_calloc(1, sizeof(Data));
    Cone * k = scs_calloc(1, sizeof(Cone));

    AMatrix * A;
	Sol sol = { 0 };
	Info info;
	static char *kwlist[] = { "shape", "Ax", "Ai", "Ap", "b", "c", "cone", "warm",
        "verbose", "normalize", "max_iters", "scale", "eps", "cg_rate", "alpha", "rho_x", NULL };
	
    /* parse the arguments and ensure they are the correct type */
#ifdef DLONG
	static char *argparse_string = "(ll)O!O!O!O!O!O!|O!O!O!lddddd";
#else
	static char *argparse_string = "(ii)O!O!O!O!O!O!|O!O!O!iddddd";
#endif
    npy_intp veclen[1];
    PyObject *x, *y, *s, *returnDict, *infoDict;
    
    d->stgs = scs_malloc(sizeof(Settings));

    /* set defaults */
	setDefaultSettings(d);

	if ( !PyArg_ParseTupleAndKeywords(args, kwargs, argparse_string, kwlist, 
        &(d->m), 
        &(d->n), 
        &PyArray_Type, &Ax,
		&PyArray_Type, &Ai, 
        &PyArray_Type, &Ap, 
        &PyArray_Type, &b, 
        &PyArray_Type, &c, 
        &PyDict_Type, &cone,
        &PyDict_Type, &warm,
        &PyBool_Type, &verbose,
        &PyBool_Type, &normalize,
        &(d->stgs->max_iters),
        &(d->stgs->scale),
        &(d->stgs->eps),
        &(d->stgs->cg_rate),
        &(d->stgs->alpha),
        &(d->stgs->rho_x)) ) {
        PySys_WriteStderr("error parsing inputs\n");
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

	/* get the typenum for the primitive scs_int and scs_float types */
	intType = getIntType();
	scs_floatType = getDoubleType();
	
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
	ps.Ax = getContiguous(Ax, scs_floatType);
	ps.Ai = getContiguous(Ai, intType);
	ps.Ap = getContiguous(Ap, intType);

    A = scs_malloc(sizeof(AMatrix));
	A->n = d->n;
    A->m = d->m;
    A->x = (scs_float *) PyArray_DATA(ps.Ax);
	A->i = (scs_int *) PyArray_DATA(ps.Ai);
	A->p = (scs_int *) PyArray_DATA(ps.Ap);
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
	ps.c = getContiguous(c, scs_floatType);
	d->c = (scs_float *) PyArray_DATA(ps.c);
	/* set b */
	if (!PyArray_ISFLOAT(b) || PyArray_NDIM(b) != 1) {
		return finishWithErr(d, k, &ps, "b must be a dense numpy array with one dimension");
	}
	if (PyArray_DIM(b,0) != d->m) {
		return finishWithErr(d, k, &ps, "b has incompatible dimension with A");
	}
	ps.b = getContiguous(b, scs_floatType);
	d->b = (scs_float *) PyArray_DATA(ps.b);
    
    if (getPosIntParam("f", &(k->f), 0, cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field f");
	}
	if (getPosIntParam("l", &(k->l), 0, cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field l");
	}
	if (getConeArrDim("q", &(k->q), &(k->qsize), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field q");
	}
	if (getConeArrDim("s", &(k->s), &(k->ssize), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field s");
	}
	if (getConeFloatArr("p", &(k->p), &(k->psize), cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field p");
	}
	if (getPosIntParam("ep", &(k->ep), 0, cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field ep");
	}
	if (getPosIntParam("ed", &(k->ed), 0, cone) < 0) {
		return finishWithErr(d, k, &ps, "failed to parse cone field ed");
	}

    d->stgs->verbose = verbose ? (scs_int) PyObject_IsTrue(verbose) : VERBOSE;
    d->stgs->normalize = normalize ? (scs_int) PyObject_IsTrue(normalize) : NORMALIZE;
    if(d->stgs->max_iters < 0) {
		return finishWithErr(d, k, &ps, "max_iters must be positive");
	}
    if(d->stgs->scale < 0) {
		return finishWithErr(d, k, &ps, "scale must be positive");
	}
    if(d->stgs->eps < 0) {
		return finishWithErr(d, k, &ps, "eps must be positive");
	}
    if(d->stgs->cg_rate < 0) {
		return finishWithErr(d, k, &ps, "cg_rate must be positive");
	}
    if(d->stgs->alpha < 0) {
		return finishWithErr(d, k, &ps, "alpha must be positive");
	}
    if(d->stgs->rho_x < 0) {
		return finishWithErr(d, k, &ps, "rho_x must be positive");
	}
	/* parse warm start if set */
    d->stgs->warm_start = WARM_START;
	if (warm) {
		d->stgs->warm_start = getWarmStart("x", &(sol.x), d->n, warm);
		d->stgs->warm_start |= getWarmStart("y", &(sol.y), d->m, warm);
		d->stgs->warm_start |= getWarmStart("s", &(sol.s), d->m, warm);
	}
	
    /* Solve! */
    scs(d, k, &sol, &info);

	/* create output (all data is *deep copied*) */
	/* x */
	/* matrix *x; */
	/* if(!(x = Matrix_New(n,1,DOUBLE))) */
	/*   return PyErr_NoMemory(); */
	/* memcpy(MAT_BUFD(x), mywork->x, n*sizeof(scs_float)); */
	veclen[0] = d->n;
	x = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.x);
    PyArray_ENABLEFLAGS((PyArrayObject *) x, NPY_ARRAY_OWNDATA);

	/* y */
	/* matrix *y; */
	/* if(!(y = Matrix_New(p,1,DOUBLE))) */
	/*   return PyErr_NoMemory(); */
	/* memcpy(MAT_BUFD(y), mywork->y, p*sizeof(scs_float)); */
	veclen[0] = d->m;
	y = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.y);
    PyArray_ENABLEFLAGS((PyArrayObject *) y, NPY_ARRAY_OWNDATA);

	/* s */
	/* matrix *s; */
	/* if(!(s = Matrix_New(m,1,DOUBLE))) */
	/*   return PyErr_NoMemory(); */
	/* memcpy(MAT_BUFD(s), mywork->s, m*sizeof(scs_float)); */
	veclen[0] = d->m;
	s = PyArray_SimpleNewFromData(1, veclen, NPY_DOUBLE, sol.s);
    PyArray_ENABLEFLAGS((PyArrayObject *) s, NPY_ARRAY_OWNDATA);

    infoDict = Py_BuildValue("{s:l,s:l,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:s}",
			"statusVal", (scs_int) info.statusVal, "iter", (scs_int) info.iter, "pobj", (scs_float) info.pobj,
			"dobj", (scs_float) info.dobj, "resPri", (scs_float) info.resPri, "resDual", (scs_float) info.resDual,
			"relGap", (scs_float) info.relGap, "resInfeas", (scs_float) info.resInfeas, "resUnbdd", (scs_float) info.resUnbdd,
			"solveTime", (scs_float) (info.solveTime), "setupTime", (scs_float) (info.setupTime),
			"status", info.status);

    returnDict = Py_BuildValue("{s:O,s:O,s:O,s:O}", "x", x, "y", y, "s", s, "info", infoDict);
	/* give up ownership to the return dictionary */
	Py_DECREF(x);
	Py_DECREF(y);
	Py_DECREF(s);
	Py_DECREF(infoDict);

	/* no longer need pointers to arrays that held primitives */
	freePyData(d, k, &ps);
	return returnDict;
}

static PyMethodDef scsMethods[] = {
    { "csolve", (PyCFunction) csolve, METH_VARARGS | METH_KEYWORDS, "Solve a convex cone problem using scs."},
    { "version", (PyCFunction)version, METH_NOARGS, "Version number for SCS."},
    { NULL, NULL, 0, NULL } /* sentinel */
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

	if (m == NULL) {
		return NULL;
    }
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
