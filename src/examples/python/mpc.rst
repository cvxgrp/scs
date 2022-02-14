.. _py_mpc:

Model Predictive Control
========================

In this example we shall demonstrate an instance of using the box :ref:`cone
<cones>`, as well as reusing a :ref:`cached workspace and using warm-starting
<warm_start>`.

In model predictive control (MPC) the control action at each time-step is
obtained by solving an optimization problem that simulates the dynamical system
over some time horizon.  Here, we consider the problem of controlling a linear,
time-invariant dynamical systems with quadratic stage costs and box constraints
on the state and action:

.. math::
  \begin{array}{ll}
    \mbox{minimize}   & x_T^T Q_T x_T + q_T^\top x_T + \sum_{t=0}^{T-1}\left( x_t^T Q x_t + u_t^T R u_t + q^\top x_t \right) \\
    \mbox{subject to} & x_{t+1} = A x_t + B u_t, \quad t = 0, \ldots T-1\\
                      & x_{\rm min} \le x_t  \le x_{\rm max}, \quad t = 0, \ldots T\\
                      & u_{\rm min} \le u_t  \le u_{\rm max}, \quad t = 0, \ldots T-1\\
                      & x_0 = \bar{x}
  \end{array}

over variables corresponding to the states :math:`x_t \in \mathbf{R}^{n}`,
:math:`t=0,\ldots,T`, and the inputs :math:`u_t \in \mathbf{R}^{m}`,
:math:`t=0,\ldots,T-1`, where :math:`A \in \mathbf{R}^{n \times n}`, :math:`B
\in  \mathbf{R}^{n \times m}` correspond to the linear dynamics and :math:`Q \in
\mathbf{R}^{n \times n}`, :math:`R \in \mathbf{R}^{m \times m}`, :math:`q \in
\mathbf{R}^{n}` correspond to the positive semidefinite quadratic cost functions
with terminal cost defined by :math:`Q_T \in  \mathbf{R}^{n \times n}` and
:math:`q_T \in \mathbf{R}^{n}`.

The problem is solved repeatedly for varying initial state :math:`\bar{x} \in
\mathbf{R}^{n}`. The upper and lower bound constraints can be expressed in
SCS using the box :ref:`cone <cones>`.

Python code to solve this is below.

.. literalinclude:: mpc.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python mpc.py > mpc.py.out

.. literalinclude:: mpc.py.out
      :language: none
