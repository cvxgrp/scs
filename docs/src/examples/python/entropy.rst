.. _py_entropy:

Maximum Entropy
===============

entropy is a well known technique for sparse linear regression.
It is obtained by adding an :math:`\ell_1` regularization term in the objective,

.. math::

  \begin{array}{ll}
    \mbox{maximize} & -\sum_i^n x_i \log x_i \\
    \mbox{subhect to} & {\bf 1}^T x_i = 1 \\
                      & Fx - g \geq 0
  \end{array}


where :math:`x \in \mathbf{R}^{n}` is the vector of parameters, :math:`A \in
\mathbf{R}^{m \times n}` is the data matrix, and :math:`\lambda > 0` is the
weighting parameter.  The problem has the following equivalent form,

.. math::

  \begin{array}{ll}
    \mbox{minimize} & -{\bf 1}^T t \\
    \mbox{subject to} & {\bf 1}^T x_i = 1 \\
                      & Fx - g \geq 0 \\
                      & \begin{bmatrix} t_i \\ x_i \\ 1 \end{bmatrix} \in K_\mathrm{exp}.
  \end{array}


In order to get a good trade-off between sparsity of the solution and quality of
the linear fit, we solve the problem for varying weighting parameter
:math:`\lambda`.  Since :math:`\lambda` enters only in the linear part of the
objective function, we can reuse the matrix factorization and enable warm
starting to reduce the computation time.

Python code to solve this is below.

.. literalinclude:: entropy.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python entropy.py > entropy.py.out

.. literalinclude:: entropy.py.out
      :language: none
