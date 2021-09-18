.. _python_interface:

Python
======

After :ref:`installing <python_install>` you can import SCS using

.. code:: python

  import scs

This module provides a single function solve with the following call signature:

.. code:: python

  sol = scs.solve(data,
                  cone,
                  use_indirect=False,
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
                  rho_x=1e-3,
                  acceleration_lookback=0,
                  acceleration_interval=1,
                  time_limit_secs=0,
                  write_data_filename=None,
                  log_csv_filename=None)

where :code:`data` is a dict containing :code:`P, A, b, c`, and :code:`cone` is
a dict that contains that :ref:`cones` information. 
The :code:`b`, and :code:`c` must be 1d numpy arrays and :code:`P` and :code:`A`
must be scipy sparse matrices in CSC format; if they are not of the proper
format, SCS will attempt to convert them. The remaining fields are explained
in :ref:`settings`.

At termination :code:`sol` is a dict with fields :code:`x, y, s, info` where
:code:`x, y, s` contains the primal-dual :ref:`solution <optimality>` or the
:ref:`certificate of infeasibility <infeasibility>`, and :code:`info` is a dict
containing the solve :ref:`info`.


Warm-starting
-------------

Warm-starting SCS with a guess of the primal-dual solution can reduce the total
solve time. This is useful, for example, when solving several similar problems
sequentially. To do this add to the :code:`data` dict passed to
:code:`scs.solve` the additional fields :code:`x`, :code:`y`, and :code:`s` (or
any subset thereof) where :code:`x` and :code:`s` correspond to numpy arrays
containing the primal solution guesses and :code:`y` corresponds to a numpy
array containing the dual solution guess.

