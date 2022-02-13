.. _py_mpc:

Model Predictive Control
========================

XXX

In this example we shall demonstrate an instance of using the box :ref:`cone
<cones>`, as well as reusing a cached workspace (matrix factorization) and using
:ref:`warm-starting <warm-starting>`.

Model 

We consider the problem of controlling a linear time-invariant dynamical system to some reference state :math:`x_r \in \mathbf{R}^{n_x}`.
To achieve this we use *constrained linear-quadratic MPC*, which solves at each time step the following finite-horizon optimal control problem

.. math::
  \begin{array}{ll}
    \mbox{minimize}   & x_T^T Q_T x_T + \sum_{t=0}^{T-1}\left( x_t^T Q x_t + u_t^T R u_t\right) \\
    \mbox{subject to} & x_{t+1} = A x_t + B u_t + d, \quad t = 0, \ldots T-1\\
                      & x_{\rm min} \le x_t  \le x_{\rm max}, \quad t = 0, \ldots T\\
                      & u_{\rm min} \le u_t  \le u_{\rm max}, \quad t = 0, \ldots T-1\\
                      & x_0 = \bar{x}
  \end{array}

over variables corresponding to the states :math:`x_t \in \mathbf{R}^{n_x}`,
:math:`t=0,\ldots,T`, and the inputs :math:`u_t \in \mathbf{R}^{n_u}`,
:math:`t=0,\ldots,T-1` where :math:`A \in \mathbf{R}^{n_x \times n_x}`,
:math:`b \in  \mathbf{r}^{n_x \times n_u}` and :math:`d \in  \mathbf{R}^{n_x}`
correspond to the linear dynamics (with a constant drift) and :math:`Q \in
\mathbf{R}^{n_x \times n_x}` :math:`R \in \mathbf{R}^{n_u \times n_u}`
correspond to the positive semidefinite quadratic cost functions with terminal
cost matrix :math:`Q_T \in  \mathbf{R}^{n_x \times n_x}`.

The problem is solved repeatedly for varying initial state :math:`\bar{x} \in
\mathbf{R}^{n_x}`. The upper and lower bound constraints can be expressed in
SCS using the box :ref:`cone <cones>`.

Python code to solve this is below.

.. literalinclude:: mpc.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python mpc.py > mpc.py.out

.. literalinclude:: mpc.py.out
      :language: none
