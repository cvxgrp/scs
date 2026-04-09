.. _python_interface:

Python
======

After :ref:`installing <python_install>` you can import SCS using

.. code:: python

  import scs

This module provides the :code:`SCS` class which is initialized using:

.. code:: python

  solver = scs.SCS(data,
                  cone,
                  use_indirect=False,
                  mkl=False,
                  apple_ldl=False,
                  gpu=False,
                  verbose=True,
                  normalize=True,
                  max_iters=int(1e5),
                  scale=0.1,
                  adaptive_scale=True,
                  eps_abs=1e-4,
                  eps_rel=1e-4,
                  eps_infeas=1e-7,
                  alpha=1.5,
                  rho_x=1e-6,
                  acceleration_lookback=10,
                  acceleration_interval=10,
                  time_limit_secs=0,
                  write_data_filename=None,
                  log_csv_filename=None)

where :code:`data` is a dict containing :code:`P, A, b, c`, and :code:`cone` is
a dict that contains the :ref:`cones` information. The :code:`cone` dict
contains keys corresponding to the cone type and values corresponding to either
the cone length or the array that defines the cone (see the third column in
:ref:`cones` for the keys and what the corresponding values represent).  The
:code:`b`, and :code:`c` entries must be 1d numpy arrays and the :code:`P` and
:code:`A` entries must be scipy sparse matrices in CSC format; if they are not
of the proper format, SCS will attempt to convert them. The
:code:`use_indirect` setting switches between the sparse direct
:ref:`linear_solver` (the default) or the sparse indirect solver. If the MKL
Pardiso direct solver for SCS is :ref:`installed <python_install>` then it can
be used by setting :code:`mkl=True`. On macOS the Apple Accelerate sparse
LDL\ :sup:`T` solver is included automatically and can be used by setting
:code:`apple_ldl=True`. If a GPU solver for SCS is :ref:`installed
<python_install>` and a GPU is available then it can be used by setting
:code:`gpu=True`. For the direct GPU solver based on cuDSS set
:code:`use_indirect=False`. The remaining fields are explained in
:ref:`settings`.

Cone dict
---------

The :code:`cone` dict supports the following keys (see :ref:`cones` for mathematical
definitions):

.. list-table::
   :header-rows: 1

   * - Key
     - Value
     - Description
   * - :code:`z`
     - :code:`int`
     - Zero cone length.
   * - :code:`l`
     - :code:`int`
     - Non-negative cone length.
   * - :code:`bu`, :code:`bl`
     - :code:`array`
     - Box cone upper/lower bounds (length :math:`\text{bsize}-1`).
   * - :code:`q`
     - :code:`list[int]`
     - Second-order cone lengths.
   * - :code:`s`
     - :code:`list[int]`
     - PSD cone matrix dimensions.
   * - :code:`cs`
     - :code:`list[int]`
     - Complex PSD cone matrix dimensions.
   * - :code:`ep`
     - :code:`int`
     - Number of primal exponential cone triples.
   * - :code:`ed`
     - :code:`int`
     - Number of dual exponential cone triples.
   * - :code:`p`
     - :code:`list[float]`
     - Power cone parameters in :math:`[-1, 1]`.

.. _python_spectral_cone_keys:

**Spectral cone keys** (require :ref:`spectral cone build <python_spectral_install>`):

.. list-table::
   :header-rows: 1

   * - Key
     - Value
     - Description
   * - :code:`d`
     - :code:`list[int]`
     - Log-determinant cone matrix dimensions.
   * - :code:`nuc_m`, :code:`nuc_n`
     - :code:`list[int]`
     - Nuclear norm cone matrix row/column dimensions (must be equal length, :math:`m_i \geq n_i`).
   * - :code:`ell1`
     - :code:`list[int]`
     - :math:`\ell_1` norm cone vector dimensions.
   * - :code:`sl_n`, :code:`sl_k`
     - :code:`list[int]`
     - Sum-of-largest-eigenvalues cone matrix dimensions and :math:`k` values (must be equal length, :math:`0 < k_i < n_i`).

Then to solve the problem call:

.. code:: python

   sol = solver.solve(warm_start=True, x=None, y=None, s=None)

where :code:`warm_start` indicates whether the solve will reuse the previous
solution as a warm-start (if this is the first solve it initializes at zero).
A good warm-start can reduce the overall number of iterations required to solve
a problem. 1d Numpy arrays :code:`x,y,s` are (optional) warm-start overrides if
you wish to set these manually rather than use solution to the last problem as
the warm-start.

At termination :code:`sol` is a dict with fields :code:`x, y, s, info` where
:code:`x, y, s` contains the primal-dual :ref:`solution <optimality>` or the
:ref:`certificate of infeasibility <infeasibility>`, and :code:`info` is a dict
containing the solve :ref:`info`.

To re-use the workspace and solve a similar problem with new :code:`b`
and / or :code:`c` data, we can update the solver using:

.. code:: python

   solver.update(b=new_b, c=new_c)  # update b and c vectors (can be None)
   solver.solve()  # solve new problem with updated b and c


