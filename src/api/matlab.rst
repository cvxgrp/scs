.. _matlab_interface:

MATLAB
======

After :ref:`installing <matlab_install>` you can call

.. code:: matlab

  [x, y, s, info] = scs(data, cones, settings)

where :code:`data` is a struct containing :code:`P, A, b, c`, :code:`P, A` must
be sparse matrices, :code:`settings` is a struct containing solver
:ref:`settings` (missing settings are set to the defaults), and :code:`cones` is
a struct that contains the :ref:`cones` information.  The :code:`cone` struct
contains members corresponding to the cone type and values corresponding to either
the cone length or the array that defines the cone (see the third column in
:ref:`cones` for the keys and what the corresponding values represent).  At
termination :code:`x, y, s` contains the primal-dual :ref:`solution
<optimality>` or the :ref:`certificate of infeasibility <infeasibility>`, and
info is a struct containing the solve :ref:`info`.

Warm-starting
-------------

Warm-starting SCS with a guess of the primal-dual solution can reduce the total
solve time. This is useful, for example, when solving several similar problems
sequentially. To do this add to the :code:`data` struct passed to :code:`scs`
the additional fields :code:`x`, :code:`y`, and :code:`s` (or any subset
thereof) where :code:`x` and :code:`s` correspond to the primal solution guesses
and :code:`y` corresponds to the dual solution guess.

Solver backends
---------------

By default SCS uses the sparse direct (LDL) solver. Alternative backends can be
selected via the :code:`settings` struct:

.. code:: matlab

  settings.use_indirect = true;  % conjugate gradient solver
  settings.dense = true;         % dense Cholesky (best for dense A)
  settings.gpu = true;           % GPU solver

Spectral cones
---------------

The MATLAB interface is compiled with :ref:`spectral cone <spectral_cones>`
support enabled by default.  This includes log-determinant, nuclear norm,
:math:`\ell_1` norm, and sum-of-largest-eigenvalues cones.  Set the
corresponding fields in the :code:`cones` struct (e.g., :code:`cones.d`,
:code:`cones.nuc_m`, :code:`cones.nuc_n`, :code:`cones.ell1`,
:code:`cones.sl_n`, :code:`cones.sl_k`) as described in the :ref:`spectral
cones documentation <spectral_cones>`.

.. note::

   Spectral cones are **experimental** — the API may change in future releases.

Workspace reuse
---------------

When solving a sequence of problems where only :code:`b` and/or :code:`c` change
(e.g., MPC, parameter sweeps), you can avoid re-factorizing by using the
workspace API:

.. code:: matlab

  % Initialize workspace (factorizes A and P)
  work = scs_init(data, cones, settings);

  % Solve
  [x, y, s, info] = scs_solve(work);

  % Update b and/or c without re-factorizing (pass [] to leave unchanged)
  scs_update(work, b_new, c_new);

  % Re-solve with the updated data
  [x, y, s, info] = scs_solve(work);

  % Warm-start a re-solve
  warm.x = x; warm.y = y; warm.s = s;
  [x, y, s, info] = scs_solve(work, warm);

  % Free workspace when done
  scs_finish(work);

This corresponds to the C API functions :code:`scs_init`, :code:`scs_solve`,
:code:`scs_update`, and :code:`scs_finish`. The workspace handle :code:`work`
records which solver backend was used, so calls to :code:`scs_solve`,
:code:`scs_update`, and :code:`scs_finish` are automatically routed to the
correct backend.
