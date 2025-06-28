.. _linear_solver:

Linear System Solver
====================

At each iteration :math:`k` SCS solves the following set of linear equations:

.. math::
  \begin{bmatrix}
  R_x + P  &  A^\top \\
  A &  -R_y   \\
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

for a particular right-hand side :math:`z^k \in \mathbf{R}^{n+m}`. The presence
of the diagonal scaling :math:`R \in \mathbf{R}^{(n+m) \times (n+m)}` matrix is
explained in :ref:`scaling`. The :math:`R_y` term is negated to make the matrix
quasi-definite; we can recover the solution to the original problem from this
modified one. Note that the matrix does not change from iteration to iteration,
which is a major advantage of this approach. Moreover, this system can be solved
approximately so long as the errors satisfy a summability condition, which
permits the use of approximate solvers that can scale to very large problems.


Available linear solvers
------------------------

Each of the below linear solvers is included in their own binary. If linking
against SCS directly, then to switch between them you must compile and link
against the right binary. If calling SCS via one of the interfaces then you can
choose between the different linear solvers using the appropriate settings.

.. _direct:

Sparse direct
^^^^^^^^^^^^^

The direct method is the default linear system solver in SCS and factorizes the
above matrix using a sparse permuted LDL factorization. Then it solves the
linear system at each iteration using the cached factors.  Since the linear
system is quasidefinite we have strong existence guarantees about the
factorization.  It relies on the external (but included) `AMD
<https://github.com/DrTimothyAldenDavis/SuiteSparse>`_ and `QDLDL
<https://github.com/oxfordcontrol/qdldl>`_ packages.

.. _mkl:

MKL Pardiso
^^^^^^^^^^^
The `Intel MKL Pardiso solver
<https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface.html>`_
is a shared-memory multiprocessing parallel direct sparse solver which is
included Intel oneAPI MKL library. This offers an alternative to the single
threaded AMD / QDLDL libraries which come bundled with SCS. Pardiso tends to be
faster than AMD / QDLDL, especially for larger problems. If MKL is installed on
your system then it is generally worth using MKL for both the blas / lapack
usage as well as the linear system solve.
Intel MKL is now available for
`free and without restrictions for everyone <https://www.intel.com/content/www/us/en/developer/articles/news/free-ipsxe-tools-and-libraries.html>`_,
though it only offers limited support for non-Intel CPUs.

.. _indirect:

Sparse indirect
^^^^^^^^^^^^^^^

The indirect method solves the above linear system approximately with a
'matrix-free' method. To do this it first reduces the system to solving

.. math::

  \begin{align}
  (R_x + P + A^\top R_y^{-1} A) x & = z^k_x + A^\top R_y^{-1} z^k_y \\
                            y & = R_y^{-1}(A x - z^k_y).
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

.. _cudss:

Sparse GPU direct method
^^^^^^^^^^^^^^^^^^^^^^^^^

This is a linear solver that uses the `cuDSS<https://developer.nvidia.com/cudss>`_
library to solve the linear system on the GPU. It is similar to the direct
solver, but uses the calls to cuDSS library to perform the analysis, numerical
(re-)factorization and subsequent solves. According to its documentation

> reordering (a major part of the analysis phase) is executed on the host,
> while symbolic factorization (another part of the analysis phase),
> numerical factorization and solve are executed on the GPU.

As the newest addition to SCS this solver is still under development and not
as well battle-tested as the other solvers.

.. _new_linear_solver:

Implementing a new linear solver
--------------------------------

In order to implement you own linear system solver, you need to implement the
struct :code:`ScsLinSysWork` that contains the workspace your solver requires,
and implement the functions in :code:`include/linsys.h` as detailed below.
See :code:`linsys` directory for examples.

.. doxygentypedef:: ScsLinSysWork

.. doxygenfile:: include/linsys.h

