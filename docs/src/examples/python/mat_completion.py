import scs
import numpy as np
import scipy as sp
from scipy import sparse

np.random.seed(1)

# The vec function as documented in api/cones
def vec(S):
    n = S.shape[0]
    S = np.copy(S)
    S *= np.sqrt(2)
    S[range(n), range(n)] /= np.sqrt(2)
    return S[np.triu_indices(n)]


# The mat function as documented in api/cones
def mat(s):
    n = int((np.sqrt(8 * len(s) + 1) - 1) / 2)
    S = np.zeros((n, n))
    S[np.triu_indices(n)] = s / np.sqrt(2)
    S = S + S.T
    S[range(n), range(n)] /= np.sqrt(2)
    return S


dim = 15  # dim x dim matrix
vlen = int(dim * (dim + 1) / 2)  # length of vector x = vec(X)

# Generate true matrix
rank = dim // 5  # low rank
X = np.random.randn(dim, rank)
X = X @ X.T

#############################################################################

# Let's first do some basic sanity checks to ensure that mat, vec are working:

# mat(vec( . )) should be identity
print(f"Should be ~ 0: {np.linalg.norm(X - mat(vec(X)))}")

# Trace( . ) should be vec(I)' vec( . )
print(f"Should be ~ 0: {np.trace(X) - vec(np.eye(dim)) @ vec(X)}")

#############################################################################

num_measurements = vlen // 2  # how many measurements are revealed

# Generate random measurement indices
measurement_idxs = np.random.choice(
    np.arange(vlen), size=num_measurements, replace=False
)

# Create A matrix
Ad = np.zeros((num_measurements, vlen))
for i in range(num_measurements):
    Ad[i, measurement_idxs[i]] = 1.0

# Noisy measurements of X
measurements = Ad @ vec(X) + 0.01 * np.random.randn(num_measurements)  # + noise

# Auxiliary data
In = sparse.eye(vlen)
Im = sparse.eye(num_measurements)
On = sparse.csc_matrix((vlen, vlen))
Onm = sparse.csc_matrix((vlen, num_measurements))

# SCS data
P = sparse.block_diag([On, sparse.eye(num_measurements)], format="csc")
A = sparse.vstack(
    [
        # zero cone
        sparse.hstack([Ad, -Im]),
        # positive semidefinite cone
        sparse.hstack([-In, Onm]),
    ],
    format="csc",
)
b = np.hstack([measurements, np.zeros(vlen)])
c = np.hstack([np.zeros(vlen + num_measurements)])

data = dict(P=P, A=A, b=b, c=c)
cone = dict(z=num_measurements, s=dim)
# Setup workspace
solver = scs.SCS(data, cone, eps_abs=1e-6, eps_rel=1e-6)
print(f"Solving for lambda = 0")
sol = solver.solve()  # lambda = 0
X_hat = mat(sol["x"][:vlen])
print(f"Error: {np.linalg.norm(X_hat - X) / np.linalg.norm(X)}")

# Solve for different values of lambda
lambdas = np.logspace(-6, 1, 11)
for lam in lambdas:
    print(f"Solving for lambda = {lam}")
    # Re-use workspace, just update the `c` vector
    c_new = np.hstack([lam * vec(np.eye(dim)), np.zeros(num_measurements)])
    solver.update(c=c_new)
    # Solve updated problem
    sol = solver.solve()  # will warm-start automatically
    X_hat = mat(sol["x"][:vlen])
    # What is the norm error?
    print(f"Error : {np.linalg.norm(X_hat - X) / np.linalg.norm(X)}")
