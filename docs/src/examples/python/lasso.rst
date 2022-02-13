.. _py_lasso:

Lasso
=====
This example demonstrates quadratic objectives, as well as
reusing a cached workspace (matrix factorization) and using warm-starting.

In the lasso the goal is to find a sparse vector that fits some measurements.
The :math:`\ell_1` norm is used as a convex surrogate for sparsity, and
a regularization parameter :math:`\lambda \geq 0` trades off sparsity and
quality of fit. Concretely the lasso is as follows

.. math::

  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} \| Ax - b \|_2^2 + \lambda \| x \|_1
  \end{array}


over variable :math:`x \in \mathbf{R}^{n}`, with data :math:`A \in
\mathbf{R}^{m \times n}` and :math:`b \in \mathbf{R}^n`.  The problem has the
following equivalent form,

.. math::

    \begin{array}{ll}
    \mbox{minimize}   & \frac{1}{2} y^T y + \lambda  {\bf 1}^T t \\
    \mbox{subject to} & y = Ax - b \\
                      & -t \le x \le t
  \end{array}

over variables :math:`x \in \mathbf{R}^{n}`, :math:`t \in \mathbf{R}^{n}`,
:math:`y \in \mathbf{R}^{m}`.  From this formulation it is straightforward to
convert it into the standard form accepted by SCS.  In order to get a good
trade-off between sparsity of the solution and quality of the linear fit, we
solve the problem for varying weighting parameter :math:`\lambda \geq 0`.  Since
:math:`\lambda` enters only in the linear part of the objective function, we can
reuse the matrix factorization and use warm starting to reduce the computation
time.

Python code to solve this is below.

.. literalinclude:: lasso.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python lasso.py > lasso.py.out

.. literalinclude:: lasso.py.out
      :language: none
