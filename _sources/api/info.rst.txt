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
     - Number of updates to the :code:`scale` parameter or to the per-row multipliers of :code:`adaptive_diag_scale` in this solve; an update of both at once is counted once. Each update requires a new factorization when a direct linear solver is used. This can be nonzero even if :code:`scale` never changes (see :ref:`updating_scale`)
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
     - Final value of the scalar :code:`scale` parameter, which can be passed as the :code:`scale` setting of a new workspace for a similar problem (see :ref:`updating_scale`). The per-row multipliers of :code:`adaptive_diag_scale` are not returned and start at one in a new workspace. When the same workspace is solved again both the :code:`scale` and the multipliers carry over from the previous solve (see :ref:`warm_start`)
   * - :code:`rejected_accel_steps`
     - :code:`scs_int`
     - Number of times an AA update was rejected by the safeguarding check (see :ref:`acceleration`)
   * - :code:`accepted_accel_steps`
     - :code:`scs_int`
     - Number of times an AA update was accepted by the safeguarding check (see :ref:`acceleration`)
   * - :code:`aa_stats`
     - :code:`AaStats`
     - Detailed AA solve diagnostics, including rejection causes, last rank, last weight norm, and last regularization
   * - :code:`lin_sys_time`
     - :code:`scs_float`
     - Total time (milliseconds) spent in the :ref:`linear system solver <linear_solver>`
   * - :code:`cone_time`
     - :code:`scs_float`
     - Total time (milliseconds) spent in the :ref:`cone projection <cones>`
   * - :code:`accel_time`
     - :code:`scs_float`
     - Total time (milliseconds) spent in the :ref:`acceleration routine <acceleration>`

Anderson acceleration statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :code:`aa_stats` field contains detailed diagnostics from the Anderson
acceleration linear solves. These counters are useful for diagnosing whether AA
is active and, when AA updates are rejected, why they were rejected.

Unlike :code:`rejected_accel_steps` and :code:`accepted_accel_steps`, which are
set to zero at the start of every solve, the :code:`aa_stats` fields belong to
the workspace and are not cleared between solves. If a workspace is reused for
several solves, the :code:`n_*` counters are totals over all the solves so far,
and the :code:`last_*` fields describe the most recent AA solve, which may have
been in an earlier solve. The exception is :code:`iter`, which counts the calls
to AA since its memory was last cleared: it restarts whenever an AA update is
rejected and whenever the :code:`scale` or the per-row multipliers are updated
(see :ref:`updating_scale`), but not at the start of a solve.

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - :code:`iter`
     - :code:`scs_int`
     - Number of calls to AA since its memory was last cleared
   * - :code:`n_accept`
     - :code:`scs_int`
     - Number of AA updates accepted by :code:`aa_apply` before safeguarding
   * - :code:`n_reject_lapack`
     - :code:`scs_int`
     - Number of AA updates rejected because the LAPACK solve failed
   * - :code:`n_reject_rank0`
     - :code:`scs_int`
     - Number of AA updates rejected because rank truncation produced rank zero
   * - :code:`n_reject_nonfinite`
     - :code:`scs_int`
     - Number of AA updates rejected because the AA weight norm was non-finite
   * - :code:`n_reject_weight_cap`
     - :code:`scs_int`
     - Number of AA updates rejected because the AA weight norm exceeded the configured cap
   * - :code:`n_safeguard_reject`
     - :code:`scs_int`
     - Number of AA updates rejected by the safeguarding check
   * - :code:`last_rank`
     - :code:`scs_int`
     - Rank used in the most recent AA solve
   * - :code:`last_aa_norm`
     - :code:`scs_float`
     - AA weight norm from the most recent AA solve, or NaN if no AA solve was attempted
   * - :code:`last_regularization`
     - :code:`scs_float`
     - Regularization value used in the most recent AA solve
