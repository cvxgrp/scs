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
   * - :code:`adaptive_diag_scale`
     - :code:`scs_int`
     - Whether to refine the metric :math:`R` per row, using the primal residual profile. Requires :code:`adaptive_scale` (silently disabled without it). See :ref:`scaling`.
     - True/False
     - 1
   * - :code:`rho_x`
     - :code:`scs_float`
     - Primal scale factor. See :ref:`scaling`.
     - :math:`(0, \infty)`
     - 1e-6 (1e-2 in single precision builds)
   * - :code:`max_iters`
     - :code:`scs_int`
     - Maximum number of iterations to run.
     - :math:`\mathbf{N}`
     - 1e5
   * - :code:`eps_abs`
     - :code:`scs_float`
     - Absolute feasibility tolerance, see :ref:`termination`. :code:`0` gives a purely relative test.
     - :math:`[0, \infty)`
     - 1e-4
   * - :code:`eps_rel`
     - :code:`scs_float`
     - Relative feasibility tolerance, see :ref:`termination`. :code:`0` gives a purely absolute test. If :code:`eps_abs` and :code:`eps_rel` are both :code:`0` the optimality test can never pass (the code tests with a strict inequality), and a problem that has a solution runs until :code:`max_iters` or :code:`time_limit_secs`.
     - :math:`[0, \infty)`
     - 1e-4
   * - :code:`eps_infeas`
     - :code:`scs_float`
     - Infeasibility tolerance (primal and dual), see :ref:`infeasibility`. :code:`0` switches the infeasibility and unboundedness tests off, so SCS never stops early with an infeasible or unbounded status.
     - :math:`[0, \infty)`
     - 1e-7
   * - :code:`alpha`
     - :code:`scs_float`
     - Douglas-Rachford relaxation parameter. See :ref:`relaxation`.
     - :math:`(0, 2)`
     - 1.5
   * - :code:`time_limit_secs`
     - :code:`scs_float`
     - Time limit for the solve in seconds (can be fractional). :code:`0` is interpreted as no limit. The clock starts when :code:`scs_solve` is called, so the setup time (e.g., the initial factorization) is not counted, although the time column of the printed output does include it. The limit is only tested when convergence is checked, every 25 iterations (:code:`CONVERGED_INTERVAL` in :code:`include/glbopts.h`), so a solve can overrun it by up to 25 iterations. When the limit is reached SCS returns its best guess with an inaccurate status, see :ref:`exit flags <exit_flags>`. It has no effect if SCS is compiled with the :code:`NO_TIMER` :ref:`compile flag <compile_flags>`. It can be changed after setup, for example to give the solve what is left of a total time budget, see :ref:`update_settings`.
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
     - 5
   * - :code:`acceleration_type_1`
     - :code:`scs_int`
     - Whether AA uses type-I (:code:`1`) or type-II (:code:`0`). See :ref:`acceleration`.
     - True/False
     - 1
   * - :code:`acceleration_regularization`
     - :code:`scs_float`
     - Tikhonov regularization :math:`\epsilon` for the AA least-squares solve. A positive value is relative: :math:`\epsilon` is the value times :math:`\|S_k\|_F \|Y_k\|_F` for type-I and times :math:`\|Y_k\|_F^2` for type-II. A negative value is absolute: :math:`\epsilon` is its magnitude, with no scaling. :code:`0` turns the regularization off. The default is tuned for type-I; type-II is typically more stable and tolerates a smaller value (e.g. :code:`1e-12`). See :ref:`acceleration`.
     - :math:`(-\infty, \infty)`
     - 1e-8
   * - :code:`acceleration_relaxation`
     - :code:`scs_float`
     - AA relaxation factor :math:`\beta`; :code:`1.0` recovers vanilla AA. See :ref:`acceleration`.
     - :math:`[0, 2]`
     - 1.0
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





.. _update_settings:

Changing settings after setup
-----------------------------

.. note::
   :code:`scs_update_settings` is new after SCS 3.3.2: SCS 3.3.2 and earlier
   do not have it, and the Python :code:`update(**settings)` below needs an
   scs-python release built on an SCS that does.

Some settings are read only while SCS iterates, so they can be changed on an
existing workspace between solves, without repeating the setup (the data
scaling and the factorization or preconditioning of the linear system). This
helps when:

* you have a **total time budget that includes the setup**.
  :code:`time_limit_secs` counts the solve only, so a solve that stops on
  the limit overruns the budget by the setup time. Set up first, measure how
  long it took, then give the solve what is left;
* you want to solve to a loose tolerance first and then tighten it, or
  change :code:`max_iters`, on the same workspace;
* you want verbose output or a CSV log for some solves and not others.

The settings that can be changed after setup are:

.. list-table::
   :header-rows: 1

   * - Setting
     - Notes
   * - :code:`time_limit_secs`
     - The clock still starts when the solve starts. :code:`0` means no limit.
   * - :code:`max_iters`
     -
   * - :code:`eps_abs`, :code:`eps_rel`, :code:`eps_infeas`
     - Also used by the :ref:`scale updates <updating_scale>`, which happen
       during the solve.
   * - :code:`alpha`
     -
   * - :code:`verbose`
     - Only the solve output: the setup output has already been printed.
   * - :code:`log_csv_filename`
     - The file is opened when each solve starts. SCS keeps its own copy of
       the name.
   * - :code:`warm_start`
     - Accepted, but in C the :code:`warm_start` argument of
       :code:`scs_solve` always replaces it.

Every other setting is used to set up the workspace and cannot be changed
after it (see the table at the end of this section). To change one of those,
set up a new workspace.

In C
^^^^

Call :code:`scs_update_settings` with a complete :ref:`ScsSettings
<ScsSettings>` struct, normally the one you passed to :code:`scs_init` with
the fields you want to change edited. The setup settings in it must still
have the values you passed to :code:`scs_init`; if any differs, or a new
value fails the checks :code:`scs_init` applies, SCS prints a message naming
the setting, changes nothing, and returns a negative value. The new values
stay in force for every later solve on the workspace until you change them
again. To spend a total budget of :code:`budget` seconds:

.. code:: c

  /* time scs_init with your own clock: the setup time SCS measures is only
   * reported in ScsInfo after a solve */
  double t0 = now_in_seconds();
  ScsWork *w = scs_init(d, k, stgs);
  double remaining = budget - (now_in_seconds() - t0);

  if (remaining > 0) {
    stgs->time_limit_secs = remaining; /* 0 would mean no limit */
    if (scs_update_settings(w, stgs) < 0) {
      /* nothing was changed; SCS printed why */
    }
    exitflag = scs_solve(w, sol, &info, 0);
  } /* else the budget is spent: do not solve */

Do not set :code:`time_limit_secs` to :code:`0` when the budget is spent:
:code:`0` means no limit. Either skip the solve or set a tiny positive value,
which stops the solve at its first convergence check, after one iteration.

In Python
^^^^^^^^^

Pass the new values to :code:`update` as keyword arguments, see
:ref:`python_interface`. Only the settings you pass change; passing a setup
setting raises an error and changes nothing.

.. code:: python

  t0 = time.perf_counter()
  solver = scs.SCS(data, cone, eps_abs=1e-6, eps_rel=1e-6)
  remaining = budget - (time.perf_counter() - t0)
  if remaining > 0:
      solver.update(time_limit_secs=remaining)
      sol = solver.solve()

How to check it worked
^^^^^^^^^^^^^^^^^^^^^^

A solve stopped by the time limit returns an inaccurate :ref:`exit flag
<exit_flags>` and its status ends in :code:`(inaccurate - reached
time_limit_secs)`. With :code:`verbose` on, :code:`scs_update_settings`
prints the values now in force, starting :code:`settings updated:`.

What it does not do
^^^^^^^^^^^^^^^^^^^

* It does not make :code:`time_limit_secs` include the setup: the limit
  still counts from the start of the solve, so subtract the setup time
  yourself as above.
* It does not make the limit exact: SCS tests the time only when it checks
  convergence, every 25 iterations, so a solve can overrun the limit by up to
  25 iterations plus the time to finish (as with any time limit). The limit
  has no effect in a build with the :code:`NO_TIMER` :ref:`compile flag
  <compile_flags>`.
* It does not change the setup time reported in the :ref:`info <info>`
  struct (:code:`scs_update` does: it reports its own time as the setup
  time).
* It does not redo any part of the setup, and it does not reset what the
  workspace carries from one solve to the next: the adapted :code:`scale`,
  the :code:`adaptive_diag_scale` multipliers and the acceleration memory
  (see :ref:`warm_start`). A solve after the call therefore starts from the
  same state it would have started from without it. Called before the first
  solve, it solves with the same settings as passing the values to
  :code:`scs_init` (the results agree up to the backend's rounding).
* It must not be called while a solve on the same workspace is running.

Settings fixed at setup, and why:

.. list-table::
   :header-rows: 1

   * - Setting
     - Why it cannot change after setup
   * - :code:`normalize`
     - The data were rescaled during setup.
   * - :code:`scale`, :code:`rho_x`
     - They are in the matrix the linear solver factorized (or
       preconditioned) during setup. :code:`scale` is compared with the value
       passed to :code:`scs_init`, not the value :code:`adaptive_scale` has
       moved it to, so passing the original struct back is accepted.
   * - :code:`adaptive_scale`, :code:`adaptive_diag_scale`
     - Setup decides from them whether the per-row multipliers of
       :code:`adaptive_diag_scale` exist (it switches
       :code:`adaptive_diag_scale` off when :code:`adaptive_scale` is off;
       the value passed to :code:`scs_init` is the one compared).
   * - :code:`acceleration_lookback`, :code:`acceleration_type_1`,
       :code:`acceleration_regularization`, :code:`acceleration_relaxation`
     - The acceleration memory was sized and configured with them during
       setup.
   * - :code:`acceleration_interval`
     - Read only by the solve, but the acceleration memory is kept from one
       solve to the next and holds steps taken at the old interval.
   * - :code:`write_data_filename`
     - The data were written during setup.
