#!/usr/bin/env python
from warnings import warn
from numpy import transpose 
from scipy import sparse
import _scs_indirect

__version__ = _scs_indirect.version()

def solve(probdata, cone, **kwargs):
    """
    solves convex cone problems
     
    @return dictionary with solution with keys:
         'x' - primal solution
         's' - primal slack solution
         'y' - dual solution
         'info' - information dictionary
    """
    if not probdata or not cone:
        raise TypeError("Missing data or cone information")

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
    if sparse.issparse(A):
        warn("Converting A to a dense matrix")
        A = A.todense()

    if sparse.issparse(b):
        b = b.todense()

    if sparse.issparse(c):
        c = c.todense()

    m, n = A.shape

    # A is stored in ROW MAJOR order, so we need to transpose:
    if kwargs.pop('gpu', False): # False by default
        import _scs_gpu
        return _scs_gpu.csolve((m, n), A.T, b, c, cone, warm, **kwargs)
    if not kwargs.pop('use_indirect', True): # True by default
        import _scs_direct
        return _scs_direct.csolve((m, n), A.T, b, c, cone, warm, **kwargs)
    return _scs_indirect.csolve((m, n), A.T, b, c, cone, warm, **kwargs)
