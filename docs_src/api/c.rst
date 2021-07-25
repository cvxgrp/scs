.. _c_interface:

C / C++
=======


.. _C_main_API:

Main solver API
---------------

The main C API is imported from the header :code:`scs.h`.

.. doxygenfunction:: scs

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

.. list-table:: Exit flags 
   :widths: 75 25
   :header-rows: 0

   * - Solved
     - 1
   * - Solved (inaccurate)
     - 2
   * - Indeterminate (numerical error) 
     - 0
   * - Unbounded
     - -1
   * - Infeasible
     - -2
   * - Indeterminate
     - -3
   * - Failed (usually data input error)
     - -4
   * - Interrupted (SIGINT)
     - -5
   * - Unbounded (inaccurate)
     - -6
   * - Infeasible (inaccurate)
     - -7
   * - Hit maximum iterations
     - -8
   * - Hit time limit
     - -9


Basic types
-----------

The most basic types are

* :code:`scs_int`: can be :code:`long` or :code:`int` if the compiler flag :code:`DLONG` is set or not
* :code:`scs_float`: can be a :code:`double` or a :code:`float` if the compiler flag :code:`SFLOAT` is set or not.


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

