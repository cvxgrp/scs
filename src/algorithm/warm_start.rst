.. _warm_start:

Caching the workspace and warm-starts
-------------------------------------

SCS supports reusing the workspace between solves so long as the data matrices
:math:`A` and :math:`P` do not change. After an initial solve the workspace can
be updated with new :math:`b` and :math:`c` vectors if desired. This can
substantially speed up subsequent solve times since we can cache important
quantities such as the matrix factorization and the data equilibration.
Moreover, SCS supports warm-starting the solver with a guess of the solution,
which can significantly reduce the total number of iterations required for
convergence.  

Re-using the workspace and warm-starting can be useful, for example, when
solving a sequence of related problems such as in :ref:`Model predictive control
<py_mpc>` or solving for the entire regularization path in the :ref:`the Lasso
<py_lasso>`.

In the :ref:`C API <c_interface>` call :code:`scs_init` once to initialize the
workspace, then use the :code:`scs_solve` in conjunction with
:code:`scs_update` to solve a sequence of problems.  Warm-starting can be done
by setting the warm-start :ref:`setting <settings>` to :code:`True` when calling
:code:`scs_solve`, and then including the guess of the solution in the :code:`x,
y, s` members of the :ref:`ScsSolution` struct, where those members correspond
to the guess of the solution in the :ref:`standard form <optimality>`.  SCS will
initialize the solver at those points and then overwrite the :ref:`ScsSolution`
struct members with the real solution at termination.

In other languages caching the workspace and warm-starting is documented in
their respective :ref:`interfaces`.

