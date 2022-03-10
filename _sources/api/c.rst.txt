.. _c_interface:

C / C++
=======


.. _C_main_API:

Main solver API
---------------

The C API is imported from the header :code:`scs.h`, available `here
<https://github.com/cvxgrp/scs/blob/master/include/scs.h>`_. This file
and :code:`scs_types.h` are the only public header files required
for downstream applications.

.. doxygenfunction:: scs_init
.. doxygenfunction:: scs_solve
.. doxygenfunction:: scs_update
.. doxygenfunction:: scs_finish

Helper functions
^^^^^^^^^^^^^^^^

This sets the :ref:`ScsSettings <ScsSettings>` struct to the default values as
specified in the :ref:`settings` page.

.. doxygenfunction:: scs_set_default_settings

If you only need to solve a single problem and not a series of related problems,
then you can call the :code:`scs` function documented here, which simply calls
the :code:`scs_init`, :code:`scs_solve`, :code:`scs_finish` sequence above.

.. doxygenfunction:: scs


..
  Lower level
  ^^^^^^^^^^^

  Under the hood the :code:`scs` function above simply calls the three functions
  below in series.  It can be useful to call :code:`scs_solve` many times for the
  same call to :code:`scs_init`. If you want to do this, for example because you
  want to cache the matrix factorization for many solves, please `contact us
  <mailto:splitting.conic.solver@gmail.com>`_, because currently that
  functionality is disabled.

  .. doxygenfunction:: scs_init

  .. doxygenfunction:: scs_solve

  .. doxygenfunction:: scs_finish


Primitive types
---------------

These are defined in header file :code:`scs_types.h`.

* :code:`scs_int`: is :code:`long` if the :ref:`compiler flag <compile_flags>` :code:`DLONG` is set, otherwise it is :code:`int`
* :code:`scs_float`: is :code:`float` if the :ref:`compiler flag <compile_flags>` :code:`SFLOAT` is set, otherwise it is :code:`double`


Input Types
-----------

The relevant input structs required by API are as follows.

.. _ScsData:

Data
^^^^

.. doxygenstruct:: ScsData
   :members:

.. _ScsMatrix:

Data Matrices
^^^^^^^^^^^^^

The matrices must be in `Compressed Sparse Column (CSC) format <https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html>`_ using zero-based indexing.
See :ref:`matrices` for more details on what SCS expects.

.. doxygenstruct:: ScsMatrix
   :members:

.. _ScsCone:

Cone
^^^^

See :ref:`cones` for more details.

.. doxygenstruct:: ScsCone
   :members:

.. _ScsSettings:

Settings
^^^^^^^^

See :ref:`settings` for details on each of these.

.. doxygenstruct:: ScsSettings
  :members:

Output Types
------------

The relevant output structs returned by SCS are as follows.

.. _ScsSolution:

Solution
^^^^^^^^

This will contain the solution as found by SCS *or* the certificate of primal or
dual infeasibility (see :ref:`termination`). If the user wants to warm-start the
solver, then the Solution struct is also used as an input to specify the
warm-start points (see :ref:`warm_start`).


.. doxygenstruct:: ScsSolution
   :members:

.. _ScsInfo:

Info
^^^^^

See :ref:`info` for details on each of these.

.. doxygenstruct:: ScsInfo
   :members:

Workspace
---------

The user should not need to interact with the :code:`ScsWork` struct,
which contains the internal workspace allocated and maintained by SCS.

