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
