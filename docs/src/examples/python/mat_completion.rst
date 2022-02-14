.. _py_mat_completion:

Low-Rank Matrix Completion
==========================

This example demonstrates how to use the :ref:`positive semidefinite cone
<sdcone>`, as well as reusing a :ref:`cached workspace and using warm-starting
<warm_start>`.

Matrix completion is the problem of filling in missing data into a partially
observed matrix where the measurements we are given have been corrupted by
Gaussian noise.  In low-rank matrix completion we have the additional prior
knowledge that the matrix we are completing is low-rank.  For simplicity we
shall also assume that the matrix we are reconstructing is symmetric positive
definite.  A famous instance of this problem is the `Netflix prize
<https://en.wikipedia.org/wiki/Netflix_Prize>`__.

Concretely, we denote by :math:`\hat X \in \mathbf{R}^{n \times n}` the true
matrix corrupted by noise, and denote by :math:`\mathcal{I}` the set of indices
(row and column pairs) from which we receive noisy observations. We shall use
the nuclear norm, denoted :math:`\| \cdot \|_*`, as a convex surrogate for rank,
which we shall trade off against the observations using regularization parameter
:math:`\lambda \geq 0`. The low-rank matrix completion problem is given by

.. math::

  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} \sum_{i,j \in \mathcal{I_n}} (X_{ij} - \hat X_{ij})^2 + \lambda \|X \|_* \\
    \mbox{subject to} & X \succeq 0
  \end{array}

over variable :math:`X \in \mathbf{R}^{n \times n}` (we use :math:`\cdot \succeq
0` to indicate membership in the symmetric positive semidefinite cone).

We can convert this into a more standard form.  First, let :math:`x =
\mathrm{vec}(X)` be the semidefinite vectorization of :math:`X` described in
:ref:`cones <sdcone>` (and concretely implemented in the code that follows).
Further, let :math:`A` be the linear operator that extracts the elements of
:math:`x` for which we have (noisy) observations, and let :math:`b = A
\mathrm{vec}(\hat X)`.  Since the nuclear norm of a positive semidefinite matrix
is given by its trace we obtain

.. math::

    \begin{array}{ll}
    \mbox{minimize}   & \frac{1}{2} y^T y + \lambda \mathrm{vec}(I_n)^\top x \\
    \mbox{subject to} & y = Ax - b \\
                      & \mathrm{mat}(x) \succeq 0
    \end{array}

over variable :math:`x \in \mathbf{R}^{n(n+1) /2}` and :math:`y \in
\mathbf{R}^{|\mathcal{I}|}`, where :math:`I_n` is the :math:`n \times n`
identity matrix.  From this formulation it is straightforward to convert it into
the standard form accepted by SCS.  The regularization parameter :math:`\lambda
\geq 0` trades off the rank of the solution and the quality of the fit, and so
we solve the problem for many choices of :math:`\lambda`.  Since :math:`\lambda`
enters only in the linear part of the objective function, we can reuse the
matrix factorization and use warm starting to reduce the computation time.

Python code to solve this is below.

.. literalinclude:: mat_completion.py
      :language: python

After following the python :ref:`install instructions <python_install>`, we can
run the code yielding output:

.. python mat_completion.py > mat_completion.py.out

.. literalinclude:: mat_completion.py.out
      :language: none
