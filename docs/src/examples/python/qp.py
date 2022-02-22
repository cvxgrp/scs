import scipy
import scs
import numpy as np

# Set up the problem data
P = scipy.sparse.csc_matrix([[3.0, -1.0], [-1.0, 2.0]])
A = scipy.sparse.csc_matrix([[-1.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
b = np.array([-1, 0.3, -0.5])
c = np.array([-1.0, -1.0])

# Populate dicts with data to pass into SCS
data = dict(P=P, A=A, b=b, c=c)
cone = dict(z=1, l=2)

# Initialize solver
solver = scs.SCS(data, cone, eps_abs=1e-9, eps_rel=1e-9)
# Solve!
sol = solver.solve()

print(f"SCS took {sol['info']['iter']} iters")
print("Optimal solution vector x*:")
print(sol["x"])

print("Optimal dual vector y*:")
print(sol["y"])
