import scs
import numpy as np
import scipy as sp
from scipy import sparse

# Generate problem data
sp.random.seed(1)

# Matrix size parameters.
n = 20
m = 10
p = 5

# Generate random problem data.
tmp = np.random.rand(n)
A = np.random.randn(m, n)
b = A.dot(tmp)
F = np.random.randn(p, n)
g = F.dot(tmp) + np.random.rand(p)


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
    sol = solver.solve()
    x = sol["x"][:n]
    # What is the norm error?
    print(f"Error : {np.linalg.norm(x_true - x) / np.linalg.norm(x_true)}")
