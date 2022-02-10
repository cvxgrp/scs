.. _py_mat_completion:

Low-Rank Matrix Completion
==========================


Matrix completion is the problem of filling in missing data into a partially
observed matrix where the measurements we are given have
been corrupted by Gaussian noise.  In low-rank matrix completion we have
the additional prior knowledge that the matrix we are completing is low-rank.
For simplicity we shall also assume that the matrix we are reconstructing is
symmetric positive definite.  A famous instance of
this problem is the `Netflix prize
<https://en.wikipedia.org/wiki/Netflix_Prize>`__.

Concretely, we denote by :math:`\hat X \in \mathcal{S}^n_+` the true matrix
corrupted by noise, and denote by :math:`\mathcal{I}` the set of indices (row
and column pairs) from which we receive noisy observations. Let :math:`X \in
\mathcal{S}^n_+` be the
variable we are solving for. The nuclear norm, denoted :math:`\| \cdot \|_*`
XXX

.. math::

  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} \sum_{i,j \in \mathcal{I_n}} (X_{ij} - \hat X_{ij})^2 + \lambda \|X \|_* \\
    \mbox{subject to} & X \succeq 0
  \end{array}

where :math:`I_n` is the :math:`n \times n` identity matrix.

We can convert this into a XX
First, let :math:`x = \mathrm{vec}(X)` be the semidefinite vectorization
of :math:`X` described in :ref:`cones` (and concretely implemented
in the code that follows). Further, let :math:`A` be the linear operator that
extracts the elements XXX and :math:`b = A \mathrm{vec}(\hat X)`.

.. math::

    \begin{array}{ll}
    \mbox{minimize}   & \frac{1}{2} y^T y + \lambda \mathrm{vec}(I)^\top x \\
    \mbox{subject to} & y = Ax - b \\
                      & \mathrm{mat}(x) \succeq 0
    \end{array}

over variable :math:`x \in \mathbf{R}^{n(n+1) /2}`.  From this formulation it is
straightforward to convert it into the standard form accepted by SCS.  In order
to get a good trade-off between rank of the solution and quality of the fit, we
solve the problem for varying weighting parameter :math:`\lambda`.  Since
:math:`\lambda` enters only in the linear part of the objective function, we can
reuse the matrix factorization and enable warm starting to reduce the
computation time.

Python code to solve this is below.

.. literalinclude:: mat_completion.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python mat_completion.py > mat_completion.py.out

.. literalinclude:: mat_completion.py.out
      :language: none
