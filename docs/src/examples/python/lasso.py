import scs
import numpy as np
import scipy as sp
from scipy import sparse

# Generate problem data
sp.random.seed(1)
np.random.seed(1)

n = 200  # Variables
m = 100  # Measurements
Ad = sparse.random(m, n, density=0.5)  # Measurement matrix
x_true = sparse.random(n, 1, density=0.1)  # True sparse vector
x_true = np.array(x_true.todense()).squeeze()

measurements = Ad @ x_true + 0.1 * np.random.randn(m)
measurements = np.array(measurements).squeeze()

# The smallest value of lambda with solution all-zeros
lambda_max = np.linalg.norm(Ad.T @ measurements, np.inf)

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
        # zero cone
        sparse.hstack([Ad, -Im, Onm.T]),
        # positive cones
        sparse.hstack([In, Onm, -In]),
        sparse.hstack([-In, Onm, -In]),
    ],
    format="csc",
)
b = np.hstack([measurements, np.zeros(n), np.zeros(n)])
c = np.zeros(2 * n + m)

data = dict(P=P, A=A, b=b, c=c)
cone = dict(z=m, l=2 * n)

print(f"Solving for lambda = 0")
# Setup workspace
solver = scs.SCS(
    data,
    cone,
    eps_abs=1e-6,
    eps_rel=1e-6,
)
sol = solver.solve()  # lambda = 0
x = sol["x"][:n]
print(f"Error : {np.linalg.norm(x_true - x) / np.linalg.norm(x_true)}")

# Solve for different values of lambda
lambdas = np.logspace(-2, np.log10(lambda_max), 11)
for lam in lambdas:
    print(f"Solving for lambda = {lam}")
    # Re-use workspace, just update the `c` vector
    c_new = np.hstack([np.zeros(n + m), lam * np.ones(n)])
    solver.update(c=c_new)
    # Solve updated problem
    sol = solver.solve()  # will warm-start automatically
    x = sol["x"][:n]
    # What is the norm error?
    print(f"Error : {np.linalg.norm(x_true - x) / np.linalg.norm(x_true)}")
