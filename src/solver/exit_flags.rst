.. _exit_flags:

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


