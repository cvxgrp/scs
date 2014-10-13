import scs 
from sys import getrefcount
from guppy import hpy

import numpy as np
from scipy import sparse
ij = np.array([[0,1,2,3],[0,1,2,3]])
A = sparse.csc_matrix(([-1.,-1.,1.,1.], ij), (4,4))
b = np.array([0.,0.,1,1])
c = np.array([1.,1.,-1,-1])
cone = {'l': 4}

print c
print b
print A
print cone
data = {'A': A, 'b':b, 'c':c}

sol = scs.solve(data, cone)
print sol

sol = scs.solve(data, cone, use_indirect = True)
print sol

sol = scs.solve(data, cone, max_iters = 500, eps = 1e-6, normalize = False)
print sol

sol = scs.solve(data, cone, max_iters = 500, eps = 1e-6, use_indirect = True)
print sol

print getrefcount(sol['x'])
h = hpy()
print h.heap()

data['x'] = sol['x']
data['y'] = sol['y']
data['s'] = sol['s']
sol = scs.solve(data, cone, max_iters = 500, eps = 1e-6)
sol = scs.solve(data, cone, max_iters = 500, eps = 1e-6, use_indirect = True)

