.. _settings:

Settings
--------

These settings control how SCS behaves during a solve.
They are set in the :ref:`ScsSettings <ScsSettings>` struct.

.. list-table::
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Permitted values
     - Default
   * - :code:`normalize`
     - :code:`scs_int`
     - Whether to perform heuristic data rescaling. See :ref:`equilibration`.
     - True/False
     - 1
   * - :code:`scale`
     - :code:`scs_float`
     - Initial dual scale factor, updated if :code:`adaptive_scale` is True. See :ref:`scaling`.
     - :math:`(0, \infty)`
     - 0.1
   * - :code:`adaptive_scale`
     - :code:`scs_int`
     - Whether to heuristically adapt dual :code:`scale` through the solve. See :ref:`scaling`.
     - True/False
     - 1
   * - :code:`rho_x`
     - :code:`scs_float`
     - Primal scale factor. See :ref:`scaling`.
     - :math:`(0, \infty)`
     - 1e-6
   * - :code:`max_iters`
     - :code:`scs_int`
     - Maximum number of iterations to run.
     - :math:`\mathbf{N}`
     - 1e5
   * - :code:`eps_abs`
     - :code:`scs_float`
     - Absolute feasibility tolerance, see :ref:`termination`.
     - :math:`(0, \infty)`
     - 1e-4
   * - :code:`eps_rel`
     - :code:`scs_float`
     - Relative feasibility tolerance, see :ref:`termination`.
     - :math:`(0, \infty)`
     - 1e-4
   * - :code:`eps_infeas`
     - :code:`scs_float`
     - Infeasibility tolerance (primal and dual), see :ref:`infeasibility`.
     - :math:`(0, \infty)`
     - 1e-7
   * - :code:`alpha`
     - :code:`scs_float`
     - Douglas-Rachford relaxation parameter. See :ref:`relaxation`.
     - :math:`(0, 2)`
     - 1.5
   * - :code:`time_limit_secs`
     - :code:`scs_float`
     - Time limit for solve run in seconds (can be fractional). :code:`0` is interpreted as no limit.
     - :math:`[0, \infty)`
     - 0
   * - :code:`verbose`
     - :code:`scs_int`
     - Whether to print solver output to stdout.
     - True/False
     - 1
   * - :code:`warm_start`
     - :code:`scs_int`
     - Set to True if you initialize the solver with a guess of the solution. See :ref:`warm_start`. This is overridden by the argument passed to :code:`scs_solve`.
     - True/False
     - 0
   * - :code:`acceleration_lookback`
     - :code:`scs_int`
     - How much memory to use for Anderson acceleration. More memory requires more time to compute but can give more reliable steps. :code:`0` disables it. See :ref:`acceleration`.
     - :math:`\mathbf{N}`
     - 10
   * - :code:`acceleration_interval`
     - :code:`scs_int`
     - Run Anderson acceleration every :code:`acceleration_interval` iterations. See :ref:`acceleration`.
     - :math:`\mathbf{N}`
     - 10
   * - :code:`write_data_filename`
     - :code:`char *`
     - If this is set the problem data is dumped to this filename.
     - Any filename
     - NULL
   * - :code:`log_csv_filename`
     - :code:`char *`
     - If this is set SCS will write csv logs of various quantities through the solver (makes the solver much slower).
     - Any filename
     - NULL



