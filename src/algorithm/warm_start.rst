.. _warm_start:

Warm-starting
-------------

SCS supports warm-starting (aka hot-starting) the solver with a guess of the
solution which can sometimes substantially improve the algorithm performance.
This is useful, for example, when solving a sequence of related problems.

In the :ref:`raw API <c_interface>` this can be done by toggling the warm-start
:ref:`setting <settings>` to :code:`True`, and then including the guess of the
solution in the :code:`x, y, s` members of the :ref:`ScsSolution` struct, where
those members correspond to the guess of the solution in the :ref:`standard form
<optimality>`.  SCS will initialize the solver at those points and then
overwrite the :ref:`ScsSolution` struct members with the real solution at
termination.

In other languages the warm-starting is documented in their respective
:ref:`interfaces`.

