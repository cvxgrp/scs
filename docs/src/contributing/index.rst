.. _contributing:

Contributing
===============
There are many ways you can contribute to SCS, and we welcome all help!
Here are some ideas (of varying difficulty):

* Answer an open `GitHub issue <https://github.com/cvxgrp/scs/issues>`_
* Improve the documentation (the website or in the code)
* Add more :ref:`examples` to the docs
* Improve the test coverage
* Improve the :ref:`Anderson acceleration <acceleration>` stability
* Improve the :ref:`data equilibration <equilibration>`
* Improve the :ref:`heuristic re-scaling <updating_scale>`
* Determine how to select the :code:`TAU_FACTOR` :ref:`term <scaling>`
* Add other new :ref:`linear system solvers <new_linear_solver>`
* Refactor the :ref:`linear solvers <linear_solver>` to only compile a single binary with all solvers
* Add :ref:`interfaces for other languages <interfaces>` (or improve the current interfaces)

If you are interested in helping out, please start by opening a `GitHub issue
<https://github.com/cvxgrp/scs/issues>`_ so we can track progress and ensure
that our priorities align.

Building the docs
-----------------
The example pages show real solver output, which is generated at build time
rather than committed, so it can never drift from the solver. Building the docs
therefore runs the C and Python examples and needs their dependencies:

.. code:: bash

  pip install sphinx sphinx-rtd-theme breathe docutils   # the docs themselves
  pip install scs numpy scipy cvxpy                      # to run the examples
  cd docs/src && make docs

A C compiler and BLAS/LAPACK are also required, for the C example. Use
:code:`make example_outputs` to regenerate just the captured output, and
:code:`PYTHON=/path/to/venv/bin/python` to point at a specific interpreter.

The one exception is :code:`examples/qp.m.out`, which is committed: refreshing
it needs a MATLAB licence, so it is updated by hand.
