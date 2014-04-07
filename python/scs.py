#!/usr/bin/env python
import _scs_direct
import _scs_indirect
from warnings import warn
from scipy import sparse


def solve(probdata, cone, opts={}, USE_INDIRECT=False):
    """
    solves convex cone problems
     
    @return dictionary with solution with keys:
         'x' - primal solution
         's' - primal slack solution
         'y' - dual solution
         'info' - information dictionary
    """
    if not 'A' in probdata or not 'b' in probdata or not 'c' in probdata:
        raise TypeError("Missing one or more of A, b, c from data dictionary")
    A = probdata['A']
    b = probdata['b']
    c = probdata['c']

    warm = {}
    if 'x' in probdata:
        warm['x'] = probdata['x']
    if 'y' in probdata:
        warm['y'] = probdata['y']
    if 's' in probdata:
        warm['s'] = probdata['s']

    if A is None or b is None or c is None:
        raise TypeError("Incomplete data specification")
    if not sparse.issparse(A):
        raise TypeError("A is required to be a sparse matrix")
    if not sparse.isspmatrix_csc(A):
        warn("Converting A to a CSC (compressed sparse column) matrix; may take a while.")
        A = A.tocsc()

    if sparse.issparse(b):
        b = b.toDense()

    if sparse.issparse(c):
        c = c.toDense()

    m, n = A.shape

    Adata, Aindices, Acolptr = A.data, A.indices, A.indptr
    if USE_INDIRECT:
        return _scs_indirect.csolve((m, n), Adata, Aindices, Acolptr, b, c, cone, opts, warm)
    else:
        return _scs_direct.csolve((m, n), Adata, Aindices, Acolptr, b, c, cone, opts, warm)
