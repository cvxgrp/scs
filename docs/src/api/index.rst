.. _api:

API
###

Loosely speaking, SCS takes data :math:`P, A, b, c, \mathcal{K}` and produces
primal-dual optimal points :math:`(x^\star, y^\star, s^\star)` or a certificate
of primal or dual :ref:`infeasibility`.

The behavior of SCS is controlled by the :ref:`settings <settings>`.  As well as
the solution SCS also returns :ref:`information <info>` about the solve process
and a status :ref:`exit flag <exit_flags>`.

**Matrices**

The matrix :math:`A` is sparse and the rows correspond to the cones of the
problem. The order of the rows of :math:`A` must be in the order that the cones
appear in the table :ref:`here <cones>`.

The matrix :math:`P` must be sparse and symmetric positive semidefinite. Only
pass the upper triangular part into SCS.


Other languages
===============
SCS is written in raw C code, with interfaces for several other languages.

:ref:`C/C++ <c_interface>`

:ref:`Python <python_interface>`

:ref:`MATLAB <matlab_interface>`

:ref:`Julia <julia_interface>`

.. toctree::
   :maxdepth: 2
   :hidden:

   c.rst
   python.rst
   matlab.rst
   julia.rst
