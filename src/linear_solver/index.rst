.. _linear_solver:

Linear System Solver
====================

At each iteration SCS must compute the following

.. math::
  z = p^k - r \tau,

where

.. math::
  \begin{align}
  p^k &= (R + M)^{-1} R \mu^k \\
  r   &= (R + M)^{-1} q
  \end{align}

the presence of the diagonal :math:`R` matrices is explained in
:ref:`scaling` (:math:`R` does *not* appear before :math:`q` in the second
expression above).  Now consider :math:`r = (R + M)^{-1} q` and recall

.. math::
  M = \begin{bmatrix}
        P  &  A^\top \\
        -A &  0   \\
      \end{bmatrix}


Specifically, we want to solve

.. math::

  \begin{bmatrix}
  r_x \\
  r_y
  \end{bmatrix}
  =
  \begin{bmatrix}
  \rho_x I + P  &  A^\top \\
  -A &  \mathrm{diag}(\rho_y)   \\
  \end{bmatrix}
  \begin{bmatrix}
  q_x \\
  q_y
  \end{bmatrix}

which is quasidefinite if we negate the bottom row:

.. math::

  \begin{bmatrix}
  r_x \\
  r_y
  \end{bmatrix}
  =
  \begin{bmatrix}
  \rho_x I + P  &  A^\top \\
  A &  -\mathrm{diag}(\rho_y)   \\
  \end{bmatrix}
  \begin{bmatrix}
  q_x \\
  -q_y
  \end{bmatrix}

A direct method factorizes the above matrix.
An indirect method can solve via:

.. math::

  \begin{align}
  (\rho_x I + P + A^\top \mathrm{diag}(\rho_y)^{-1} A) r_x & = q_x - A^\top \mathrm{diag}(\rho_y)^{-1} q_y \\
                            r_y & = \mathrm{diag}(\rho_y)^{-1}(A z_x + q_y).
  \end{align}


Available linear solvers
------------------------

Each of the below linear solvers is included in their own binary. If linking
against SCS directly, then to switch between them you must compile and link
against the right binary. If using SCS via one of the interfaces then you can
choose between the different linear solvers using the appropriate settings.

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

The indirect method solves the above linear system approximately with a
'matrix-free' method. To do this it first reduces the system as described above
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
however for very large problems the indirect solver can be faster.

.. _gpu_indirect:

GPU indirect method
^^^^^^^^^^^^^^^^^^^

The above linear solvers all run on CPU. We also have support for a GPU version
of the indirect solver, where the multiplies are all performed on the GPU.

.. _new_lin_solver:

Implementing a new linear solver
--------------------------------

In order to implement you own linear system solver, you need to implement the
struct :code:`ScsLinSysWork` that contains the workspace your solver requires,
and implement the functions in :code:`include/linsys.h`:

.. doxygenfile:: include/linsys.h

