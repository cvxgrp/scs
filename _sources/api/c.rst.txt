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

Exit flags
-----------
The integer values that SCS can return are documented below. 

.. list-table::
   :widths: 50 10 40
   :header-rows: 1

   * - Status
     - Value
     - SCS constant name
   * - Solved
     - 1
     - :code:`SCS_SOLVED`
   * - Solved (inaccurate)
     - 2
     - :code:`SCS_SOLVED_INACCURATE`
   * - Unfinished (should never be returned) 
     - 0
     - :code:`SCS_UNFINISHED`
   * - Unbounded
     - -1
     - :code:`SCS_SOLVED`
   * - Infeasible
     - -2
     - :code:`SCS_INFEASIBLE`
   * - Indeterminate
     - -3
     - :code:`SCS_INDETERMINATE`
   * - Failed (typically data input error)
     - -4
     - :code:`SCS_FAILED`
   * - Interrupted (SIGINT)
     - -5
     - :code:`SCS_SIGINT`
   * - Unbounded (inaccurate)
     - -6
     - :code:`SCS_UNBOUNDED`
   * - Infeasible (inaccurate)
     - -7
     - :code:`SCS_INFEASIBLE`
   * - Hit maximum iterations
     - -8
     - :code:`SCS_MAX_ITERS`
   * - Hit time limit
     - -9
     - :code:`SCS_TIME_LIMIT`

Basic types
-----------

The most basic types are

* :code:`scs_int`: can be :code:`long` or :code:`int` if the compiler flag :code:`DLONG` is set or not
* :code:`scs_float`: can be a :code:`double` or a :code:`float` if the compiler flag :code:`SFLOAT` is set or not.


Types
-----------

The relevant structures used in the API are

Data
^^^^

.. doxygenstruct:: ScsData
   :members:

Matrices
^^^^^^^^

The matrices are defined in `Compressed Sparse Column (CSC) format <https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html>`_ using zero-based indexing.

.. doxygenstruct:: ScsMatrix
   :members:

Cone
^^^^

.. doxygenstruct:: ScsCone
   :members:


Settings
^^^^^^^^

.. doxygenstruct:: ScsSettings
  :members:

Solution
^^^^^^^^

.. doxygenstruct:: ScsSolution
   :members:

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

