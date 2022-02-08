.. _py_lasso:

Lasso
=====

Lasso is a well known technique for sparse linear regression.
It is obtained by adding an :math:`\ell_1` regularization term in the objective,

.. math::
     \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} \| Ax - b \|_2^2 + \lambda \| x \|_1
  \end{array}


where :math:`x \in \mathbf{R}^{n}` is the vector of parameters, :math:`A \in
\mathbf{R}^{m \times n}` is the data matrix, and :math:`\lambda > 0` is the
weighting parameter.  The problem has the following equivalent form,

.. math::
     \begin{array}{ll}
    \mbox{minimize}   & \frac{1}{2} y^T y + \lambda \boldsymbol{1}^T t \\
    \mbox{subject to} & y = Ax - b \\
                      & -t \le x \le t
  \end{array}


In order to get a good trade-off between sparsity of the solution and quality of
the linear fit, we solve the problem for varying weighting parameter
:math:`\lambda`.  Since :math:`\lambda` enters only in the linear part of the
objective function, we can reuse the matrix factorization and enable warm
starting to reduce the computation time.

Python code to solve this is below.

.. literalinclude:: lasso.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python lasso.py > lasso.py.out

.. literalinclude:: lasso.py.out
      :language: none
