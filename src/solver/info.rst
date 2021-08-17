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
   * - :code:`status_val`
     - :code:`scs_int`
     - Status integer :ref:`exit flag <exit_flags>`
   * - :code:`scale_updates`
     - :code:`scs_int` 
     - Number of updates to the scale parameter
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
   * - :code:`res_unbbd_a`
     - :code:`scs_float` 
     - Dual infeasibility residual involving :math:`A` (see :ref:`termination conditions <termination>`)
   * - :code:`res_unbdd_p`
     - :code:`scs_float` 
     - Dual infeasibility residual involving :math:`P` (see :ref:`termination conditions <termination>`)
   * - :code:`setup_time`
     - :code:`scs_float` 
     - Time taken for setup (milliseconds) 
   * - :code:`solve_time`
     - :code:`scs_float` 
     - Time taken for solve (milliseconds) 
   * - :code:`scale`
     - :code:`scs_float` 
     - Final scale parameter (useful for initializing next solve)

