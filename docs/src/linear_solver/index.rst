.. _linear_solver:

Linear System Solver
====================

At each iteration :math:`k` SCS solves the following set of linear equations:

.. math::
  \begin{bmatrix}
  \rho_x I + P  &  A^\top \\
  A &  -\mathrm{diag}(\rho_y)   \\
  \end{bmatrix}
  \begin{bmatrix}
  x \\
  y
  \end{bmatrix}
  =
  \begin{bmatrix}
  z^k_x \\
  z^k_y
  \end{bmatrix}

for a particular right-hand side :math:`z^k \in \mathbf{R}^{n+m}`
(the presence of the diagonal scaling :math:`\rho` terms is explained in
:ref:`scaling`). Note that the matrix does not change from iteration to
iteration, which is a major advantage of this approach. Moreover, this system
can be solved approximately so long as the errors satisfy a summability
condition, which permits the use of approximate solvers that can scale to
very large problems.



Available linear solvers
------------------------

Each of the below linear solvers is included in their own binary. If linking
against SCS directly, then to switch between them you must compile and link
against the right binary. If using SCS via one of the interfaces then you can
choose between the different linear solvers using the appropriate settings.

.. _direct:

Sparse direct method
^^^^^^^^^^^^^^^^^^^^

The direct method is the default linear system solver in SCS and factorizes the
above matrix using a sparse permuted LDL factorization. Then it solves the
linear system at each iteration using the cached factors.  Since the linear
system is quasidefinite we have strong existence guarantees about the
factorization.  It relies on the external `AMD
<https://github.com/DrTimothyAldenDavis/SuiteSparse>`_ and `QDLDL
<https://github.com/oxfordcontrol/qdldl>`_ packages.

.. _indirect:

Sparse indirect method
^^^^^^^^^^^^^^^^^^^^^^

The indirect method solves the above linear system approximately with a
'matrix-free' method. To do this it first reduces the system to solving

.. math::

  \begin{align}
  (\rho_x I + P + A^\top \mathrm{diag}(\rho_y)^{-1} A) r_x & = q_x - A^\top \mathrm{diag}(\rho_y)^{-1} q_y \\
                            r_y & = \mathrm{diag}(\rho_y)^{-1}(A z_x + q_y).
  \end{align}

then solves the positive definite system using using `conjugate gradients
<https://en.wikipedia.org/wiki/Conjugate_gradient_method>`_.  Each iteration of
CG requires one multiply each of sparse matrices :math:`P, A, A^\top`.  The
system is solved up to some tolerance, which is tuned to ensure that the overall
algorithm converges. The tolerance decays with iteration :math:`k` like
:math:`O(1/k^\gamma)` where :math:`\gamma > 1` and is determined by the constant
:code:`CG_RATE` (defaults to :math:`1.5`).

The indirect method has the advantage of not requiring an expensive
factorization but typically is slower on a per-iteration basis. In most cases
the factorization is relatively cheap so the direct method is the default,
however for very large problems the indirect solver can be faster.

.. _gpu_indirect:

Sparse GPU indirect method
^^^^^^^^^^^^^^^^^^^^^^^^^^

The above linear solvers all run on CPU. We also have support for a GPU version
of the indirect solver, where the matrix multiplies are all performed on the
GPU.

.. _new_linear_solver:

Implementing a new linear solver
--------------------------------

In order to implement you own linear system solver, you need to implement the
struct :code:`ScsLinSysWork` that contains the workspace your solver requires,
and implement the functions in :code:`include/linsys.h` as detailed below.
See :code:`linsys` directory for examples.

.. doxygentypedef:: ScsLinSysWork

.. doxygenfile:: include/linsys.h

