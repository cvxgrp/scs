.. _py_entropy:

Maximum Entropy
===============

This example demonstrates an instance of using the exponential :ref:`cone
<cones>`.  In this problem we want find the maximum entropy point inside a
convex polytope, ie, to solve

.. math::

  \begin{array}{ll}
    \mbox{maximize} & -\sum_i^n x_i \log x_i \\
    \mbox{subhect to} & {\bf 1}^T x = 1 \\
                      & Ax - b \geq 0
  \end{array}


over variable :math:`x \in \mathbf{R}^{n}`, where :math:`A \in
\mathbf{R}^{m \times n}` and :math:`b \in \mathbf{R}^m` are data.  The problem
has the following equivalent form,

.. math::

  \begin{array}{ll}
    \mbox{minimize} & -{\bf 1}^T t \\
    \mbox{subject to} & {\bf 1}^T x = 1 \\
                      & Ax - b \geq 0 \\
                      & \begin{bmatrix} t_i \\ x_i \\ 1 \end{bmatrix} \in \mathcal{K}_\mathrm{exp}, \quad i=1,\ldots,n,
  \end{array}

over variables :math:`x \in \mathbf{R}^{n}`, :math:`t \in \mathbf{R}^{n}` and
where :math:`\mathcal{K}_\mathrm{exp} \subset \mathbf{R}^3` denotes the
exponential cone. 

Python code to solve this is below.

.. literalinclude:: entropy.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python entropy.py > entropy.py.out

.. literalinclude:: entropy.py.out
      :language: none
