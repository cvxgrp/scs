import scs
import numpy as np
import scipy as sp
from scipy import sparse

np.random.seed(1)


class MPC(object):
    """Model Predictive Contoller using SCS."""

    def __init__(self, Ad, Bd, Q, R, q, QT, qT, xmin, xmax, umin, umax, T):
        # State and action dimension
        self.nx, self.nu = Bd.shape

        # Stack variables as follows:
        # [x_0, x_1, ..., x_{T-1}, x_T, u_0, u_1, ..., u_{T-1}]

        # Quadratic objective
        P = sparse.block_diag(
            [sparse.kron(sparse.eye(T), Q), QT, sparse.kron(sparse.eye(T), R)],
            format="csc",
        )
        # Linear objective
        c = np.hstack([np.kron(np.ones(T), q), qT, np.zeros(T * self.nu)])
        # Linear dynamics
        Ax = sparse.kron(sparse.eye(T + 1), -sparse.eye(self.nx)) + sparse.kron(
            sparse.eye(T + 1, k=-1), Ad
        )
        Bu = sparse.kron(
            sparse.vstack([sparse.csc_matrix((1, T)), sparse.eye(T)]), Bd
        )
        Aeq = sparse.hstack([Ax, Bu])

        # Will update this later with initial state
        beq = np.zeros((T + 1) * self.nx)

        # Box constraints on state and action
        Aineq = sparse.eye((T + 1) * self.nx + T * self.nu)

        box_lower = np.hstack(
            [np.kron(np.ones(T + 1), xmin), np.kron(np.ones(T), umin)]
        )
        box_upper = np.hstack(
            [np.kron(np.ones(T + 1), xmax), np.kron(np.ones(T), umax)]
        )

        A = sparse.vstack(
            [
                # zero cone
                Aeq,
                # Box cone {(t, s) | -t l <= s <= t u }
                sparse.csc_matrix((1, (T + 1) * self.nx + T * self.nu)),
                -Aineq,
            ],
            format="csc",
        )

        # Box cone add constraint t=1
        self.b = np.hstack([beq, 1, np.zeros((T + 1) * self.nx + T * self.nu)])

        data = dict(P=P, A=A, b=self.b, c=c)
        cone = dict(z=(T + 1) * self.nx, bu=box_upper, bl=box_lower)
        # Create an SCS object
        self.solver = scs.SCS(data, cone, eps_abs=1e-5, eps_rel=1e-5)

    def control(self, x0):
        # Overwrite b with new initial state
        self.b[: self.nx] = -x0
        # Update b
        self.solver.update(b=self.b)
        sol = self.solver.solve()  # will warm-start automatically
        if sol["info"]["status"] != "solved":
            raise ValueError("SCS failed to solve the problem.")

        # Return first action
        return sol["x"][-T * self.nu : -(T - 1) * self.nu]


# States dimension
nx = 20
# Control dimension
nu = 5

# State dynamics matrices
Ad = 0.1 * np.random.randn(nx, nx)  # State -> State
Bd = np.random.randn(nx, nu)  # Control -> State

# Cost matrices
Q = sparse.eye(nx)  # State
QT = 10 * Q  # Terminal State
R = 0.1 * sparse.eye(nu)  # Control

# Linear cost vector
q = 0.1 * np.random.randn(nx)
qT = q

# Initial state
x0 = 10 * np.random.randn(nx)

# Prediction horizon
T = 30

# Bounds on state
xmax = np.inf * np.ones(nx)
xmin = -np.inf * np.ones(nx)

# Bounds on control
umax = np.ones(nu)
umin = -np.ones(nu)

# Initialize Model Predictive Controller
mpc = MPC(Ad, Bd, Q, R, q, QT, qT, xmin, xmax, umin, umax, T)

# Simulate in closed loop
nsteps = 10  # Number of steps
for i in range(nsteps):
    # Get control action
    u = mpc.control(x0)
    print(f"Control action: {u}")

    # Apply first control input and update to next state
    x0 = Ad @ x0 + Bd @ u + 0.01 * np.random.normal(nx, 1)  # + noise
    x0 = np.maximum(np.minimum(x0, xmax), xmin)  # Bound to xmin, xmax
