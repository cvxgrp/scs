.. _exit_flags:

Exit flags
-----------
The integer values that SCS can return are documented below and are defined
in the :code:`include/scs.h` file.

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
   * - Indeterminate: reached :code:`max_iters` or :code:`time_limit_secs`
       at an iterate that gives no estimate of a solution or of a
       certificate; x, y and s are NaN (see :ref:`below <exit_flags_limit>`)
     - -3
     - :code:`SCS_INDETERMINATE`
   * - Failed: a data input error, an internal error such as a failed
       linear solve (the printed message says which; an internal error asks
       you to open an issue), or an iterate that is not finite (NaN or Inf)
       at :code:`max_iters` or :code:`time_limit_secs`
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


.. _exit_flags_limit:

At max_iters or time_limit_secs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When SCS stops at :code:`max_iters` or :code:`time_limit_secs` before its
:ref:`termination criteria <termination>` hold, it reports what the last
iterate supports. SCS works with a homogeneous iterate :math:`(x, y, s, \tau,
\kappa)`: a solution estimate is :math:`(x, y, s) / \tau`, which needs
:math:`\tau > 0`, and a certificate estimate is :math:`y / (-b^\top y)` or
:math:`(x, s) / (-c^\top x)`, which needs :math:`b^\top y < 0` or
:math:`c^\top x < 0`. The tests below are exact comparisons, taken in this
order, so the status depends only on the last iterate:

.. list-table::
   :widths: 40 15 45
   :header-rows: 1

   * - Last iterate
     - Exit flag
     - Returned
   * - Not finite (a NaN or an infinite entry, for example after an
       overflow)
     - -4 (:code:`SCS_FAILED`)
     - x, y, s and the info fields listed for -3 NaN; status
       :code:`failure`
   * - :math:`\kappa > \tau` and :math:`b^\top y < 0` or :math:`c^\top x
       < 0`
     - -7 or -6
     - The infeasibility certificate :math:`y / (-b^\top y)` (x and s NaN),
       or the unboundedness certificate :math:`(x, s) / (-c^\top x)` (y NaN)
   * - :math:`\tau > 0`
     - 2
     - The solution estimate :math:`(x, y, s) / \tau` and its objective and
       residuals
   * - Otherwise: :math:`\tau = 0` and :math:`b^\top y \geq 0`,
       :math:`c^\top x \geq 0`
     - -3 (:code:`SCS_INDETERMINATE`)
     - x, y, s NaN, and :code:`pobj`, :code:`dobj`, :code:`gap`,
       :code:`res_pri`, :code:`res_dual`, :code:`res_infeas`,
       :code:`res_unbdd_a`, :code:`res_unbdd_p` and :code:`comp_slack` NaN

Each status except :code:`failure` ends in :code:`(inaccurate - reached
max_iters)` or :code:`(inaccurate - reached time_limit_secs)`, for example
:code:`indeterminate (inaccurate - reached max_iters)`.

The last iterate is never built from an accelerated step that the
:ref:`acceleration <acceleration>` safeguard rejected: AA makes no step at the
last iteration, and the time limit is tested only at a convergence check,
which reads only iterates the safeguard has accepted.

The estimates returned with 2, -6 and -7 do not meet the tolerances. A
certificate estimate rests on the sign of :math:`b^\top y` or :math:`c^\top
x` at the last iterate alone, which on a run whose iterates diverged can
change from one iteration to the next.

SCS returns :code:`SCS_INDETERMINATE` when the cone projection has set
:math:`\tau` to zero and neither :math:`b^\top y` nor :math:`c^\top x` is
negative. This happens in the first few iterations and on runs whose iterates
diverged. As for an inaccurate status, more iterations or more time can help
(see :ref:`help_indeterminate`). A warm start from the NaN solution starts
from zero, as a cold start does (see :ref:`warm_start`).
