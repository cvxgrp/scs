:orphan:

API
===

SCS takes data :math:`P, A, b, c, \mathcal{K}` and produces primal-dual
:ref:`optimal <optimality>` points :math:`(x^\star, y^\star, s^\star)` or a
certificate of primal or dual :ref:`infeasibility`. The pages below describe
the supported cones, the matrix input format, the settings, the information
returned about a solve and the exit flags, followed by the C and Python
interfaces. The other language interfaces are documented on the website at
https://www.cvxgrp.org/scs/api/.

.. toctree::
   :maxdepth: 2

   /api/cones
   /api/matrices
   /api/settings
   /api/info
   /api/exit_flags
   /api/compile_flags
   /api/c
   /api/python
