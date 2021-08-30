.. _linear_solver:

Linear System Solver
====================

At each iteration SCS solves a system

.. math::

  \begin{bmatrix}
  r_x \\
  r_y
  \end{bmatrix}
  =
  \begin{bmatrix}
  R_x + P  &  A^\top \\
  A &  -R_y   \\
  \end{bmatrix}
  \begin{bmatrix}
  q_x \\
  -q_y
  \end{bmatrix}

(the presence of the diagonal :math:`R` matrices is explained in
:ref:`scaling`).

Available linear solvers
------------------------

.. _direct:

Direct method
^^^^^^^^^^^^^

The direct method is the default linear system solver in SCS and factorizes the
above matrix using a sparse permuted LDL factorization. Then it solves the
linear system at each iteration using the cached factors.  It relies on the
external `AMD <https://github.com/DrTimothyAldenDavis/SuiteSparse>`_ and `QDLDL
<https://github.com/oxfordcontrol/qdldl>`_ packages.

.. _indirect:

Indirect method
^^^^^^^^^^^^^^^

TODO(fix setting)

The indirect method can be enabled via the :code:`use_indirect` :ref:`setting
<settings>` and solves the above linear system approximately with a
'matrix-free' method. To do this it first reduces the system to

.. math::

  \begin{align}
  (R_x + P + A^\top R_y^{-1} A) r_x & = q_x - A^\top R_y^{-1} q_y \\
                            r_y & = R_y^{-1}(A z_x + q_y).
  \end{align}

then solves the positive definite system using using `conjugate gradients
<https://en.wikipedia.org/wiki/Conjugate_gradient_method>`_.  Each iteration of
CG requires one multiply each of :math:`P, A, A^\top`.  The system is solved up
to some tolerance, which is tuned to ensure that the overall algorithm
converges. The tolerance decays with iteration :math:`k` like
:math:`O(1/k^\gamma)` where :math:`\gamma > 1` and is determined by the constant
:code:`CG_RATE` (defaults to :math:`1.5`).

The indirect method has the advantage of not requiring an expensive
factorization but typically is slower on a per-iteration basis. In most cases
the factorization is relatively cheap so the direct method is the default,
however for very large problems it can be faster.

.. _gpu_indirect:

GPU indirect method
^^^^^^^^^^^^^^^^^^^

TODO

The above linear solvers all run on CPU. We also have support for a GPU version
of the indirect solver, where the multiplies are all performed on the GPU.

.. _new_lin_solver:

Implementing a new linear solver
--------------------------------
TODO
