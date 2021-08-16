.. _optimality:

Optimality conditions
---------------------

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


.. _infeasiblity:

Certificate of infeasibility
----------------------------

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
--------------------

Optimality
^^^^^^^^^^
The iterates produced by SCS *always* satisfy the conic constraints :math:`s \in
\mathcal{K}, y \in \mathcal{K}^*, s \perp y`.  Therefore to say that a problem
is solved we need to check if the primal residual, dual residual, and duality
gap are all below a certain tolerance. Specifically, SCS terminates when it has
found :math:`x \in \mathbf{R}^n`, :math:`s \in \mathbf{R}^m`, and :math:`y \in
\mathbf{R}^m` that satisfy the following residual bounds.

Primal residual:

.. math::
  r_p := \|Ax + s - b\|_\infty \leq \epsilon_\mathrm{abs} + \epsilon_\mathrm{rel} \max(\|Ax\|_\infty, \|s\|_\infty, \|b\|_\infty)

Dual residual:

.. math::
  r_d := \|Px + A^\top y - c \|_\infty \leq \epsilon_\mathrm{abs} + \epsilon_\mathrm{rel} \max(\|Px\|_\infty, \|A^\top y\|_\infty, \|c\|_\infty)

Duality gap:

.. math::
  r_g := |x^\top Px + c^\top x + b^\top y| \leq \epsilon_\mathrm{abs} + \epsilon_\mathrm{rel} \max(|x^\top P x|, |c^\top x|, |b^\top y|),

where :math:`\epsilon_\mathrm{abs}>0` and :math:`\epsilon_\mathrm{rel}>0` are user defined
:ref:`settings <settings>`.

.. _infeasibility:

Infeasibility
^^^^^^^^^^^^^

Since the conic constraints are always guaranteed by the iterates (i.e., :math:`s \in
\mathcal{K}, y \in \mathcal{K}^*, s \perp y`), SCS
declares a problem primal infeasible when it finds :math:`y \in \mathbf{R}^m` that satisfies

Primal infeasibility residual:

.. math::
  b^\top y = -1, \quad \|A^\top y\|_\infty < \epsilon_\mathrm{infeas}.

Similarly, SCS declares dual infeasibility when it finds :math:`x \in
\mathbf{R}^n`, :math:`s \in \mathbf{R}^m` that satisfy

Dual infeasibility residual:

.. math::
  c^\top x = -1, \quad  \max(\|P x\|_\infty, \|A x + s\|_\infty) < \epsilon_\mathrm{infeas}

where :math:`\epsilon_\mathrm{infeas}` is a user-defined :ref:`setting <settings>`.


