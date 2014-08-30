from __future__ import print_function
import scs
import numpy as np
from numpy import random
from scipy import sparse
# import ecos # use to verify result

#############################################
#  Uses scs to solve:                       #
#                                           # 
#        min. ||x||_1                       # 
#        s.t. Ax = b                        #
#                                           #     
#############################################

i = 50
p = 50 * i
q = 25 * i
random.seed(0)

A = sparse.rand(q, p, 0.01)
Ae = sparse.hstack([A, sparse.csc_matrix((q, p))], format="csc")
h = np.zeros(2 * p)
b = random.randn(q)
bt = np.hstack([b, h])  # in cone formulation
c = np.hstack([np.zeros(p), np.ones(p)])
I = sparse.eye(p)
G = sparse.vstack([sparse.hstack([I, -I]), sparse.hstack([-I, -I])], format("csc"))
At = sparse.vstack([Ae, G], format="csc")  # in cone formulation

# setup cone problem data
data = {'A': At, 'b': bt, 'c': c}
m = np.shape(data['A'])[0]
n = np.shape(data['A'])[1]
cone = {'l': 2 * p, 'f': q}
opts = {'normalize': True}

# indirect solver:
sol_i = scs.solve(data, cone, use_indirect=True, **opts)
# direct solver:
sol_d = scs.solve(data, cone, **opts)

# use ECOS to verify result (if installed)
# dims = {'l': 2 * p}
# sol_e = ecos.solve(c, G, h, dims, Ae, b)
# # calculate relative error:
# print(np.linalg.norm(sol_i['x'] - sol_e['x'])/np.linalg.norm(sol_e['x']))

# perturb input data, test warm-starting:
data['b'] = np.hstack([b + 0.005 * random.randn(q), h])
# put warm start guesses in data dictionary:
data['x'] = sol_i['x']
data['y'] = sol_i['y']
data['s'] = sol_i['s']
# indirect solver, warm start:
opts['cg_rate'] = 2
sol_i_warm = scs.solve(data, cone, use_indirect=True, **opts)
# direct solver, warm start:
sol_d_warm = scs.solve(data, cone, **opts)
