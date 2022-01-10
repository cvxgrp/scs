.. _exit_flags:

Exit flags
-----------
The integer values that SCS can return are documented below and are defined
in the :code:`'include/scs.h` file.

.. list-table::
   :widths: 50 10 40
   :header-rows: 1

   * - Status
     - Value
     - SCS constant name
   * - Solved to desired tolerance
     - 1
     - :code:`SCS_SOLVED`
   * - Did not reach desired accuracy, returning best guess of solution
     - 2
     - :code:`SCS_SOLVED_INACCURATE`
   * - Unfinished (never returned, only used internally)
     - 0
     - :code:`SCS_UNFINISHED`
   * - Primal unbounded / Dual infeasible (to desired tolerance)
     - -1
     - :code:`SCS_UNBOUNDED`
   * - Primal infeasible / Dual unbounded (to desired tolerance)
     - -2
     - :code:`SCS_INFEASIBLE`
   * - Indeterminate (numerical errors when recovering solution) DEPRECATED
     - -3
     - :code:`SCS_INDETERMINATE`
   * - Failed (usually a data input error)
     - -4
     - :code:`SCS_FAILED`
   * - Interrupted (received SIGINT)
     - -5
     - :code:`SCS_SIGINT`
   * - Did not reach desired accuracy, returning best guess of certificate of primal unboundedness
     - -6
     - :code:`SCS_UNBOUNDED_INACCURATE`
   * - Did not reach desired accuracy, returning best guess of certificate of primal infeasibility
     - -7
     - :code:`SCS_INFEASIBLE_INACCURATE`

