.. _info:

Return information
------------------
When SCS terminates it will return an :ref:`ScsInfo <scsinfo>` struct containing
the following fields.


.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - :code:`iter`
     - :code:`scs_int`
     - Number of iterations taken
   * - :code:`status`
     - :code:`char *`
     - Status string (e.g., 'solved')
   * - :code:`lin_sys_solver`
     - :code:`char *`
     - Linear system solver used
   * - :code:`status_val`
     - :code:`scs_int`
     - Status integer :ref:`exit flag <exit_flags>`
   * - :code:`scale_updates`
     - :code:`scs_int`
     - Number of updates to the scale parameter  (see :ref:`updating_scale`)
   * - :code:`pobj`
     - :code:`scs_float`
     - Primal objective
   * - :code:`dobj`
     - :code:`scs_float`
     - Dual objective
   * - :code:`res_pri`
     - :code:`scs_float`
     - Primal residual (see :ref:`termination conditions <termination>`)
   * - :code:`res_dual`
     - :code:`scs_float`
     - Dual residual (see :ref:`termination conditions <termination>`)
   * - :code:`gap`
     - :code:`scs_float`
     - Absolute duality gap  (see :ref:`termination conditions <termination>`)
   * - :code:`res_infeas`
     - :code:`scs_float`
     - Primal infeasibility residual (see :ref:`termination conditions <termination>`)
   * - :code:`res_unbdd_a`
     - :code:`scs_float`
     - Dual infeasibility residual involving :math:`A` (see :ref:`termination conditions <termination>`)
   * - :code:`res_unbdd_p`
     - :code:`scs_float`
     - Dual infeasibility residual involving :math:`P` (see :ref:`termination conditions <termination>`)
   * - :code:`comp_slack`
     - :code:`scs_float`
     - Complementary slackness (:math:`s^\top y`), should be very close to zero
   * - :code:`setup_time`
     - :code:`scs_float`
     - Time taken for setup (milliseconds)
   * - :code:`solve_time`
     - :code:`scs_float`
     - Time taken for solve (milliseconds)
   * - :code:`scale`
     - :code:`scs_float`
     - Final scale parameter, useful for initializing next solve (see :ref:`updating_scale`)
   * - :code:`rejected_accel_steps`
     - :code:`scs_int`
     - Number of times an AA update was rejected by the safeguarding check (see :ref:`acceleration`)
   * - :code:`accepted_accel_steps`
     - :code:`scs_int`
     - Number of times an AA update was accepted by the safeguarding check (see :ref:`acceleration`)
   * - :code:`lin_sys_time`
     - :code:`scs_float`
     - Total time (milliseconds) spent in the :ref:`linear system solver <linear_solver>`
   * - :code:`cone_time`
     - :code:`scs_float`
     - Total time (milliseconds) spent in the :ref:`cone projection <cones>`
   * - :code:`accel_time`
     - :code:`scs_float`
     - Total time (milliseconds) spent in the :ref:`aceleration routine <acceleration>`

