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

Basic types
-----------

The most basic types are

* :code:`scs_int`: is :code:`long` if the :ref:`compiler flag <compile_flags>` :code:`DLONG` is set, otherwise it is :code:`int`
* :code:`scs_float`: is :code:`float` if the :ref:`compiler flag <compile_flags>` :code:`SFLOAT` is set, otherwise it is :code:`double`


Types
-----------

The relevant structures used in the API are

.. _ScsData:

Data
^^^^

.. doxygenstruct:: ScsData
   :members:

.. _ScsMatrix:

Matrices
^^^^^^^^

The matrices are defined in `Compressed Sparse Column (CSC) format <https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html>`_ using zero-based indexing.

.. doxygenstruct:: ScsMatrix
   :members:

.. _ScsCone:

Cone
^^^^

.. doxygenstruct:: ScsCone
   :members:

.. _ScsSettings:

Settings
^^^^^^^^

.. doxygenstruct:: ScsSettings
  :members:

.. _ScsSolution:

Solution
^^^^^^^^

.. doxygenstruct:: ScsSolution
   :members:

.. _ScsInfo:

Info
^^^^^

.. doxygenstruct:: ScsInfo
   :members:

Work
^^^^^

The user should not need to interact with the :code:`ScsWork` struct,
which contains the internal workspace allocated and maintained by SCS.

.. doxygenstruct:: ScsWork
   :members:

