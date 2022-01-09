import scipy
import scs
import numpy as np

# Set up the problem data
P = scipy.sparse.csc_matrix([[3., -1.], [-1., 2.]])
A = scipy.sparse.csc_matrix([[-1., 1.], [1., 0.], [0., 1.]])
b = np.array([-1, 0.3, -0.5])
c = np.array([-1., -1.])

# Populate dicts with data to pass into SCS
data = dict(P=P, A=A, b=b, c=c)
cone = dict(l=len(b));

# Solve!
sol = scs.solve(data, cone, eps_abs=1e-9, eps_rel=1e-9)

print(f"SCS took {sol['info']['iter']} iters")
print("Optimal solution vector x*:");
print(sol['x'])

print("Optimal dual vector y*:");
print(sol['y'])
