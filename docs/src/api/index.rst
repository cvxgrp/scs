.. _api:

API
===

.. toctree::
   :maxdepth: 2
   :hidden:

   cones.rst
   matrices.rst
   settings.rst
   info.rst
   exit_flags.rst
   compile_flags.rst

Loosely speaking, SCS takes data :math:`P, A, b, c, \mathcal{K}` and produces
primal-dual :ref:`optimal <optimality>` points :math:`(x^\star, y^\star,
s^\star)` or a certificate of primal or dual :ref:`infeasibility`.  The
supported cones are documented :ref:`here <cones>`. The input format for the
data matrices is documented :ref:`here <matrices>`.  The behavior of SCS is
controlled by the :ref:`settings <settings>`.  As well as the :ref:`solution
<ScsSolution>`, SCS also returns :ref:`information <info>` about the solve
process and a status :ref:`exit flag <exit_flags>`.

.. _interfaces:

Interfaces
----------

SCS is written in raw C code, with interfaces for several other languages.

:ref:`C/C++ <c_interface>`

:ref:`Python <python_interface>`

:ref:`MATLAB <matlab_interface>`

:ref:`Julia <julia_interface>`

:ref:`R <r_interface>`

:ref:`Ruby <ruby_interface>`

:ref:`JavaScript / WebAssembly <javascript_interface>`

.. toctree::
   :maxdepth: 2
   :hidden:

   c.rst
   python.rst
   matlab.rst
   julia.rst
   r.rst
   ruby.rst
   javascript.rst
