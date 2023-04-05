import scs
import numpy as np
from scipy import sparse

# Generate problem data
np.random.seed(1)

# Matrix size parameters
n = 50  # Number of variables
p = 20  # Number of constraints

# Generate random problem data
tmp = np.random.rand(n)
tmp /= np.sum(tmp)
Ad = np.random.randn(p, n)
bd = 0.5 * Ad.dot(tmp) + 0.01 * np.random.rand(p)

# Build the A, b rows corresponding to the exponential cone
A_exp = sparse.lil_matrix((3 * n, 2 * n))
b_exp = np.zeros(3 * n)
for i in range(n):
    A_exp[i * 3, i] = -1  # t
    A_exp[i * 3 + 1, i + n] = -1  # x
    b_exp[i * 3 + 2] = 1

A = sparse.vstack(
    [
        # zero cone
        sparse.hstack([sparse.csc_matrix((1, n)), np.ones((1, n))]),
        # positive cone
        sparse.hstack([sparse.csc_matrix((p, n)), -Ad]),
        # exponential cones
        A_exp,
    ],
    format="csc",
)
b = np.hstack([1, -bd, b_exp])
c = np.hstack([-np.ones(n), np.zeros(n)])

# SCS data
data = dict(A=A, b=b, c=c)
# ep is exponential cone (primal), with n triples
cone = dict(z=1, l=p, ep=n)

# Setup workspace
solver = scs.SCS(
    data,
    cone,
)
sol = solver.solve()
x_scs = sol["x"][-n:]

# Verify solution with CVXPY
try:
    import cvxpy as cp
except ModuleNotFoundError:
    print("This example requires CVXPY installed to run.")
    raise

x = cp.Variable(shape=n)
obj = cp.Maximize(cp.sum(cp.entr(x)))
constraints = [cp.sum(x) == 1, Ad @ x >= bd]
prob = cp.Problem(obj, constraints)
prob.solve(solver=cp.ECOS)
x_cvxpy = x.value

print(f"CVXPY optimal value is:", prob.value)
print(f"Solution norm difference: {np.linalg.norm(x_scs - x_cvxpy, np.inf)}")
