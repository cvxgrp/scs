.. _solver:

Solver
===============

SCS is derived as Douglas-Rachford splitting applied to a homogeneous embedding
of the quadratic cone program. The high level algorithm is as follows,
from an initial :math:`w^0` for :math:`k=0,1,\ldots` do

.. math::
  \begin{align}
  \tilde u^{k+1} &= (I + \mathcal{Q})^{-1} w^k\\
  u^{k+1} &= \Pi_{\mathcal{C}_+}(2\tilde u^{k+1} - w^k)\\
  w^{k+1} &= w^k + u^{k+1} - \tilde u^{k+1}
  \end{align}

and :math:`u^k \rightarrow u^\star` from which we can recover the optimal solution
or the certificate of infeasibility. The algorithm consists of three steps.
The first step involves :ref:`solving a linear system` of equations:

.. math::
  \begin{bmatrix} I + P & A^T  \\ A & -I \end{bmatrix}\begin{bmatrix} x \\ y
  \end{bmatrix} = \begin{bmatrix} \mu_x \\ -\mu_y \end{bmatrix},

The second step is the Euclidean projection onto a convex :ref:`cone <cones>`, ie,

.. math::
  \begin{array}{ll}
  \mbox{minimize} & \| z - z_0\|_2^2 \\ 
  \mbox{subject to } & z \in \mathcal{K}
  \end{array}

Most cone projection operators have relatively simple projection operators.

:ref:`Termination Criteria <termination>`

:ref:`Cones <cones>`

:ref:`Settings <settings>`

Acceleration
------------

Warm-starting
-------------

SCS supports warm-starting the solver with a guess of the solution.

.. toctree::
   :maxdepth: 2
   :hidden:

   cones.rst
   termination.rst
   settings.rst
   compile_flags.rst
   info.rst
   scale.rst

