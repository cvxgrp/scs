.. _help:

Help
====

Currently the easiest way to get support is to file a `GitHub issue
<https://github.com/cvxgrp/scs/issues>`_ or to `email us
<mailto:splitting.conic.solver@gmail.com>`_ (which may be slower).

If you have a problem that SCS struggles to solve you can set the
:code:`write_data_filename` field in the :ref:`settings <settings>` and SCS will
dump a file containing the problem data to disk under that filename. Zip the
file and `email it to us <mailto:splitting.conic.solver@gmail.com>`_ or attach
it to a GitHub issue. This makes it much easier for us to reproduce the problem.

A common cause of issues is not linking :ref:`BLAS/LAPACK libraries
<blas_lapack>` correctly. If you are having this issue please search for
resources on installing and linking these libraries first. You can try `OpenBLAS
<https://www.openblas.net/>`_ if you need a BLAS library.

.. _troubleshooting:

Troubleshooting
---------------

SCS returns ``solved_inaccurate``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This means SCS ran out of iterations (or hit the time limit) before reaching
the desired accuracy. The returned solution may still be usable, but the
residuals exceed the requested tolerances. Try:

* **Increase** :code:`max_iters` (default 100,000). Some problems simply need
  more iterations to converge.
* **Loosen tolerances**: increase :code:`eps_abs` and :code:`eps_rel` if you
  don't need high accuracy.
* **Enable normalization**: set :code:`normalize=1` (the default). This
  equilibrates the problem data and generally improves convergence. If you
  turned it off, try turning it back on.
* **Try adaptive scaling**: set :code:`adaptive_scale=1` (the default). This
  automatically tunes the internal :code:`scale` parameter during the solve.

SCS returns ``infeasible`` or ``unbounded`` unexpectedly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you believe your problem is feasible and bounded but SCS disagrees:

* **Check cone ordering**: the rows of :math:`A` and :math:`b` **must** follow
  the cone order listed in :ref:`cones`. This is the most common source of
  incorrect results.
* **Check your data**: verify that :math:`A`, :math:`P`, :math:`b`, :math:`c`
  are correct. A single wrong sign or transposed matrix can make a feasible
  problem appear infeasible.
* **Tighten the infeasibility tolerance**: decrease :code:`eps_infeas` (default
  ``1e-7``). SCS may be declaring infeasibility prematurely on a
  nearly-infeasible problem.
* **Dump the problem data**: set :code:`write_data_filename` to save the exact
  data SCS sees, then inspect it.

SCS is very slow
^^^^^^^^^^^^^^^^

* **Try the indirect solver** for large, sparse problems. The direct solver
  factorizes the KKT matrix, which can be expensive for large problems. The
  indirect (conjugate gradient) solver avoids this factorization. See
  :ref:`indirect`.
* **Try MKL Pardiso** if you have Intel MKL installed. It offers a parallel
  direct solver that can be significantly faster than the default QDLDL
  for large problems. See :ref:`mkl`.
* **Check problem scaling**: if the data has entries that vary by many orders
  of magnitude, SCS will struggle. Try to rescale your problem so that the
  entries of :math:`A`, :math:`b`, and :math:`c` are roughly of order 1.
  Setting :code:`normalize=1` helps but cannot fix extreme scaling.
* **Reduce logging overhead**: if :code:`log_csv_filename` is set, SCS
  computes full residuals at every iteration, which is expensive. Remove it
  unless you need per-iteration logging.

Numerical warnings or NaN in output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Poorly conditioned data** is the most common cause. Check whether your
  constraint matrix :math:`A` or objective matrix :math:`P` is near-singular
  or has very large condition number.
* **Try adjusting** :code:`scale`: the default is ``0.1``. For some problems
  a larger value (e.g., ``1.0`` or ``10.0``) or smaller value (e.g.,
  ``0.01``) works better. Setting :code:`adaptive_scale=1` lets SCS tune
  this automatically.
* **Ensure** :math:`P` **is symmetric positive semidefinite**: SCS requires
  that only the upper triangle of :math:`P` is provided. If :math:`P` has
  negative eigenvalues the problem is non-convex and SCS cannot solve it.

Complementary slackness warning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The warning ``large complementary slackness residual`` means that
:math:`s^\top y` is not close to zero at the solution. This is usually
harmless and indicates a near-degenerate or poorly-scaled problem. The primal
and dual solutions are still valid if the other residuals are small.

