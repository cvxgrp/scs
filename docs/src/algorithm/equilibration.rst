.. _equilibration:

Data equilibration
==================

Data scale has a large impact on the convergence of first-order algorithms like
SCS in practice. With that in mind SCS implements a heuristic data equilibration
(or normalization) scheme that attempts to whiten the data which is enabled by
the :code:`normalize` :ref:`setting <settings>`.

At a high level the equilibration scheme (which is performed once, right at
the start of the solve) produces diagonal matrices :math:`E \in
\mathbf{R}^n`, :math:`D \in \mathbf{R}^m` and scalar :math:`\sigma > 0` that rescale
the data in-place, so that the tuple :math:`(P, A, b, c)` becomes
:math:`(\hat P, \hat A, \hat b, \hat c)` where

.. math::
  \hat P = EPE,\quad  \hat A  = DAE,\quad  \hat c  = \sigma Ec,\quad  \hat b = \sigma Db

Essentially:

.. math::

  \begin{bmatrix}
  \hat P & \hat A^\top & \hat c\\
  \hat A & 0 & \hat b \\
  \hat c^\top & \hat b^\top & 0
  \end{bmatrix}
  =
  \begin{bmatrix}
  E & 0 & 0 \\
  0 & D & 0 \\
  0 & 0 & \sigma
  \end{bmatrix}
  \begin{bmatrix}
  P & A^\top & c\\
  A & 0 & b \\
  c^\top & b^\top & 0
  \end{bmatrix}
  \begin{bmatrix}
  E & 0 & 0 \\
  0 & D & 0 \\
  0 & 0 & \sigma
  \end{bmatrix}
  =
  \begin{bmatrix}
  EPE & EA^\top D & \sigma Ec\\
  DAE & 0 & \sigma Db \\
  \sigma c^\top E & \sigma b^\top D & 0
  \end{bmatrix}


In other words :math:`D` rescales the rows of :math:`A, b` and :math:`E`
rescales the columns of :math:`A` and the rows / cols of :math:`P, c`.
Equilibrating matrices is a classic problem in linear algebra, and SCS simply
implements :code:`NUM_RUIZ_PASSES` (default 25) steps of Ruiz equilibration
followed by :code:`NUM_L2_PASSES` (default 1) steps of :math:`\ell_2`
equilibration on the matrix

.. math::
  \begin{bmatrix}
  P & A^\top & c\\
  A & 0 & b \\
  c^\top & b^\top & 0
  \end{bmatrix}


The main complication of the equilibration is that :math:`D` has to maintain
cone memberships, ie, it needs to satisfy for each sub-cone :math:`\mathcal{K}`

.. math::
   s \in \mathcal{K} \Leftrightarrow D_{\mathcal{K}} s \in \mathcal{K}

To do this we restrict :math:`D_\mathcal{K} = d_\mathcal{K} I_\mathcal{K}`, ie,
:math:`D` is constant diagonal within each cone. We take the scalar value to be
the :math:`\ell_\infty` norm of the calculated in-cone values for Ruiz
equilibration and the average :math:`\ell_2` norm of the in-cone values for
:math:`\ell_2` equilibration.

Originally the code supported separate :code:`primal_scale` and
:code:`dual_scale` parameters scaling :math:`c` and :math:`b`, but the latest
version of SCS has set both of these to the same value of :math:`\sigma`.

Solution
--------

This equilibration changes the fixed point of SCS, but we can recover the
solutions to the original problem from the equilibrated solution. Denote by
:math:`(\hat x, \hat y, \hat s)` the solution to the equilibrated
problem, then the solution to the original problem is given by

.. math::
   x = E\hat x / \sigma, \quad y = D \hat y / \sigma, \quad s = D^{-1} \hat s / \sigma

with objective terms

.. math::
  c^\top x =  \hat c^\top \hat x / \sigma^2, \quad
  b^\top y =  \hat b^\top \hat y / \sigma^2, \quad
  x^\top P x =  \hat x^\top \hat P \hat x / \sigma^2

The same modifications hold for the certificates of infeasibility.

Residuals
---------

The equilibration also changes the residuals and when running SCS we are
interested in detecting convergence using the original problem residuals, rather
than the equilibration problem residuals. Using the above results we can compute
the residuals of the original problem as follows (under any norm).

Primal residual:

.. math::

  r_p = \|A x + s - b\| = (1/\sigma) \| D^{-1} (\hat A \hat x + \hat s + \hat b)\|

Dual residual:

.. math::

  r_d = \|P x + A^\top y + c\| = (1/\sigma) \|E^{-1} (\hat P \hat x + \hat A^\top \hat y + \hat c) \|

Duality gap:

.. math::

  r_g = |x^\top P x + b^\top y + c^\top x|  = (1/\sigma^2) |\hat x^\top \hat P \hat x + \hat b^\top \hat y + \hat c^\top \hat x|


Infeasibility certificates
--------------------------
Similar relationships hold for the infeasibility certificate terms:

.. math::

  \|Ax + s\| = (1/\sigma) \|D^{-1} (\hat A \hat x + \hat s) \|,
  \quad
  \|Px\| = (1/\sigma) \|E^{-1} \hat P \hat x \|,
  \quad
  \|A^\top y \| = (1/\sigma) \|E^{-1} \hat A^\top \hat y \|
