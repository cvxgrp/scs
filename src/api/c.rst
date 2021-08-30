.. _c_interface:

C / C++
=======


.. _C_main_API:

Main solver API
---------------

The main C API is imported from the header :code:`scs.h`.

.. doxygenfunction:: scs

Helper functions
^^^^^^^^^^^^^^^^

.. doxygenfunction:: scs_set_default_settings

..
  lower level:
  .. doxygenfunction:: scs_init

  |

  .. doxygenfunction:: scs_solve

  |

  .. doxygenfunction:: scs_finish

Primitive types
---------------

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
See :ref:`matrices` for more details.

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

.. doxygenstruct:: ScsWork
   :members:

