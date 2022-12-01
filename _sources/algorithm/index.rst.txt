.. _algorithm:

Algorithm
===============

.. toctree::
   :maxdepth: 2
   :hidden:

   equilibration.rst
   scale.rst
   warm_start.rst
   acceleration.rst
   relaxation.rst


SCS applies Douglas-Rachford splitting to a homogeneous embedding
of the quadratic cone program. The high level algorithm is as follows,
from an initial :math:`w^0` for :math:`k=0,1,\ldots` do

.. math::
  \begin{align}
  \tilde u^{k+1} &= (I + \mathcal{Q})^{-1} w^k\\
  u^{k+1} &= \Pi_{\mathcal{C}_+}(2\tilde u^{k+1} - w^k)\\
  w^{k+1} &= w^k + u^{k+1} - \tilde u^{k+1}
  \end{align}

which yields :math:`w^k \rightarrow u^\star + \mathcal{Q}(u^\star)` where
:math:`0 \in \mathcal{Q}(u^\star) + N_{\mathcal{C}_+}(u^\star)`, and :math:`u^k
\rightarrow u^\star` from which we can recover the :ref:`optimal solution
<optimality>` or a :ref:`certificate of infeasibility <infeasibility>` (such a
:math:`u^\star` is guaranteed to exist). Note that this is slightly different to
the formulation in the original `paper
<https://web.stanford.edu/~boyd/papers/scs.html>`_, and corresponds to the
latest algorithm described `here <https://arxiv.org/abs/2004.02177>`_.  The
algorithm consists of three steps.  The first step involves :ref:`solving a
linear system <linear_solver>` of equations:

.. math::
  \begin{bmatrix} I + P & A^T  \\ A & -I \end{bmatrix}\begin{bmatrix} x \\ y
  \end{bmatrix} = \begin{bmatrix} \mu_x \\ -\mu_y \end{bmatrix}

The second step is the Euclidean projection onto a convex :ref:`cone <cones>`, ie,

.. math::
  \begin{array}{ll}
  \mbox{minimize} & \| z - z_0\|_2^2 \\
  \mbox{subject to } & z \in \mathcal{K}
  \end{array}

over variable :math:`z`. Many cones of interest have relatively simple
projection operators.

.. _optimality:

Optimality conditions
+++++++++++++++++++++
SCS solves problems of the form:

.. math::
  \begin{array}{lcr}
  \begin{array}{ll}
  \mbox{minimize} & (1/2)x^\top  P x + c^\top  x\\
  \mbox{subject to} &  Ax + s = b\\
    & s \in \mathcal{K}
  \end{array}
  &\quad&
  \begin{array}{ll}
  \mbox{maximize} & -(1/2)x^\top  P x - b^\top  y\\
  \mbox{subject to} &  Px + A^\top y + c = 0\\
    & y \in \mathcal{K}^*
  \end{array}
  \end{array}



When strong duality holds, the Karush-Kuhn-Tucker (KKT) conditions
are `necessary and sufficient for optimality <https://web.stanford.edu/~boyd/cvxbook/>`_. They are given by

.. math::
  Ax + s = b, \quad Px + A^\top y + c = 0,\quad s \in \mathcal{K}, \quad y \in \mathcal{K}^*,\quad  s \perp y.

These are primal feasibility, dual feasibility, primal and dual cone membership,
and complementary slackness.  The complementary slackness condition is
equivalent to a zero *duality gap* condition at any optimal point, that is
for :math:`(x,y,s)` that satisfy the KKT conditions we have

.. math::
  s\perp y \ \Leftrightarrow \ c^\top x + b^\top y + x^\top P x = 0.


.. _infeasibility:

Certificate of infeasibility
++++++++++++++++++++++++++++

On the other hand, if no optimal solution exists then SCS is able to robustly detect primal or dual infeasibility.
Any :math:`y \in \mathbf{R}^m` that satisfies

.. math::
  A^\top y = 0,\  y \in \mathcal{K}^*,\ b^\top y < 0

acts a certificate that the quadratic cone program is primal infeasible (dual unbounded). On the other hand, if we can find :math:`x \in \mathbf{R}^n` such
that

.. math::
  Px = 0,\ -Ax\in \mathcal{K},\ c^\top x < 0

then this is a certificate that the problem is dual infeasible (primal
unbounded).


.. _termination:

Termination criteria
++++++++++++++++++++

Optimality
----------

The iterates produced by SCS *always* satisfy the conic constraints :math:`s \in
\mathcal{K}, y \in \mathcal{K}^*, s \perp y`.  Therefore to say that a problem
is solved we need to check if the primal residual, dual residual, and duality
gap are all below a certain tolerance. Specifically, SCS terminates when it has
found :math:`x \in \mathbf{R}^n`, :math:`s \in \mathbf{R}^m`, and :math:`y \in
\mathbf{R}^m` that satisfy the following residual bounds.

* Primal residual:

.. math::
  r_p := \|Ax + s - b\|_\infty \leq \epsilon_\mathrm{abs} + \epsilon_\mathrm{rel} \max(\|Ax\|_\infty, \|s\|_\infty, \|b\|_\infty)

* Dual residual:

.. math::
  r_d := \|Px + A^\top y + c \|_\infty \leq \epsilon_\mathrm{abs} + \epsilon_\mathrm{rel} \max(\|Px\|_\infty, \|A^\top y\|_\infty, \|c\|_\infty)

* Duality gap:

.. math::
  r_g := |x^\top Px + c^\top x + b^\top y| \leq \epsilon_\mathrm{abs} + \epsilon_\mathrm{rel} \max(|x^\top P x|, |c^\top x|, |b^\top y|),

where :math:`\epsilon_\mathrm{abs}>0` and :math:`\epsilon_\mathrm{rel}>0` are
user defined :ref:`settings <settings>`. The :math:`\ell_\infty` norm
here can be changed to other norms by changing the definition of :code:`NORM` in
the :code:`include/glbopts.h` file.

.. _infeasibility_cert:

Infeasibility
-------------

Since the conic constraints are always guaranteed by the iterates (:math:`s \in
\mathcal{K}, y \in \mathcal{K}^*, s \perp y`), SCS
declares a problem **primal infeasible (dual unbounded)** when it finds :math:`y \in
\mathbf{R}^m` that satisfies

.. math::
  b^\top y = -1, \quad \|A^\top y\|_\infty < \epsilon_\mathrm{infeas}.

Similarly, SCS declares a problem **dual infeasible (primal unbounded)** when it finds
:math:`x \in \mathbf{R}^n`, :math:`s \in \mathbf{R}^m` that satisfy

.. math::
  c^\top x = -1, \quad  \max(\|P x\|_\infty, \|A x + s\|_\infty) < \epsilon_\mathrm{infeas}

where :math:`\epsilon_\mathrm{infeas} > 0` is a user-defined :ref:`setting
<settings>`.  The :math:`\ell_\infty` norm here can be changed to other norms by
changing the definition of :code:`NORM` in the :code:`include/glbopts.h` file.

In some rare cases a problem is both primal and dual infeasible. In this case
SCS will return one of the two above certificates, whichever one it finds
first. However, in that case the interpretation of infeasibility in one space
being equivalent to unboundedness in the dual space does not hold.

