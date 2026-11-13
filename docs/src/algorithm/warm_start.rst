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

The workspace also keeps the :code:`scale` parameter and the per-row
multipliers of :code:`adaptive_diag_scale` as the previous solve left them (see
:ref:`Dynamic scale updating <updating_scale>`), together with the linear
system set up for them (the factorization, for a direct solver). When
:code:`adaptive_scale` is enabled, a later solve on the same workspace, with or
without a warm-start and with or without a call to :code:`scs_update` in
between, therefore starts from these adapted values and not from the
:code:`scale` that was passed to :code:`scs_init`. Only the bookkeeping of the
update procedure starts again in each solve: the samples it averages and the
count of iterations since the last update. So the first update in a solve is
made no earlier than iteration 100, and :code:`scale_updates` in the
:ref:`info <info>` struct counts the updates of that solve only. One
consequence is that solving the same problem twice on one workspace need not
take the same number of iterations both times. To start again from the
original :code:`scale`, create a new workspace with :code:`scs_init`.

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

