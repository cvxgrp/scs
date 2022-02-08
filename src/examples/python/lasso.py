import scs
import numpy as np
import scipy as sp
from scipy import sparse

# Generate problem data
sp.random.seed(1)
n = 500
m = 5000
Ad = sparse.random(m, n, density=0.5)
x_true = np.multiply(
    (np.random.rand(n) > 0.8).astype(float), np.random.randn(n)
) / np.sqrt(n)
b = Ad.dot(x_true) + 0.5 * np.random.randn(m)
lambdas = np.linspace(1, 10, 21)

# Auxiliary data
In = sparse.eye(n)
Im = sparse.eye(m)
On = sparse.csc_matrix((n, n))
Onm = sparse.csc_matrix((n, m))

# SCS data
P = sparse.block_diag([On, sparse.eye(m), On], format="csc")
q = np.zeros(2 * n + m)
A = sparse.vstack(
    [
        sparse.hstack([Ad, -Im, Onm.T]),
        sparse.hstack([In, Onm, -In]),
        sparse.hstack([-In, Onm, -In]),
    ],
    format="csc",
)
l = np.hstack([b, np.zeros(n), np.zeros(n)])

c = np.zeros(2 * n + m)
data = dict(P=P, A=A, b=l, c=c)
cone = dict(z=m, l=2 * n)
# Setup workspace
solver = scs.SCS(
    data,
    cone,
)
sol = solver.solve()  # lambda = 0

# Solve problem for different values of lambda parameter
# Update linear cost
for lam in lambdas:
    print(f"Solving for lambda = {lam}")
    c_new = np.hstack([np.zeros(n + m), lam * np.ones(n)])
    solver.update(c=c_new)
    sol = solver.solve()
