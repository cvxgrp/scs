.. _acceleration:

Acceleration
============

SCS includes Anderson acceleration (AA), which can be used to speed up
convergence. AA is a quasi-Newton method for the acceleration of fixed point
iterations and can dramatically speed up convergence in practice, especially if
higher accuracy solutions are desired. However, it can also cause severe
instability of the solver and so should be used with caution. It is an open
research question how to best implement AA in practice to ensure good
performance across all problems and we welcome any :ref:`contributions
<contributing>` in that direction!


Mathematical details
--------------------

The discussion here is taken from section 2 of our `paper
<https://web.stanford.edu/~boyd/papers/nonexp_global_aa1.html>`__.
Consider the problem of finding a fixed point of the function :math:`f:
\mathbf{R}^n \rightarrow \mathbf{R}^n`, i.e.,

.. math::
  \text{Find } x \in \mathbf{R}^n \text{ such that } x = f(x).

In our case :math:`f` corresponds to one step of :ref:`Douglas-Rachford
splitting <algorithm>` and the iterate is the :math:`w^k` vector, which
converges to a fixed point of the DR operator. At a high level AA, from initial
point :math:`x_0` and max memory :math:`m` (corresponding to the
:code:`acceleration_lookback` setting), works as follows:

.. math::
  \begin{array}{l}
  \text{For } k=0, 1, \dots \\
  \quad \text{Set } m_k=\min\{m, k\} \\
  \quad \text{Select weights } \alpha_j^k \text{ based on the last } m_k \text{
    iterations satisfying } \sum_{j=0}^{m_k}\alpha_j^k=1 \\
  \quad \text{Set } x^{k+1}=\sum_{j=0}^{m_k}\alpha_j^kf(x^{k-m_k+j})
  \end{array}

In other words, AA produces an iterate that is the linear combination of the
last :math:`m_k + 1` outputs of the map.  Thus, the main challenge is in
choosing the weights :math:`\alpha \in \mathbf{R}^{m_k+1}`. There are two ways
to choose them, corresponding to type-I and type-II AA (named for the type-I and
type-II Broyden updates). We shall present type-II first.

Type-II AA
""""""""""

Define the residual :math:`g: \mathbf{R}^n \rightarrow \mathbf{R}^n` of
:math:`f` to be :math:`g(x) = x - f(x)`. Note that any fixed point
:math:`x^\star` satisfies :math:`g(x^\star) = 0`.
In type-II AA the weights are selected by solving a small least squares problem.

.. math::
  \begin{array}{ll}
  \mbox{minimize} & \|\sum_{j=0}^{m_k}\alpha_j g(x^{k-m_k+j})\|_2^2\\
  \mbox{subject to} & \sum_{j=0}^{m_k}\alpha_j=1,
  \end{array}

More explicitly, we can reformulate the above as follows:

.. math::
  \begin{array}{ll}
  \mbox{minimize} & \|g_k-Y_k\gamma\|_2,
  \end{array}

with variable :math:`\gamma=(\gamma_0,\dots,\gamma_{m_k-1}) \in
\mathbf{R}^{m_k}`. Here :math:`g_i=g(x^i)`,
:math:`Y_k=[y_{k-m_k}~\dots~y_{k-1}]` with :math:`y_i=g_{i+1}-g_i` for each
:math:`i`, and :math:`\alpha` and :math:`\gamma` are related by
:math:`\alpha_0=\gamma_0`, :math:`\alpha_i=\gamma_i-\gamma_{i-1}` for
:math:`1\leq i\leq m_k-1` and :math:`\alpha_{m_k}=1-\gamma_{m_k-1}`.

Assuming that :math:`Y_k` is full column rank, the solution
:math:`\gamma^k` to the above is given by :math:`\gamma^k=(Y_k^\top
Y_k)^{-1}Y_k^\top g_k`, and hence by the relation between :math:`\alpha^k` and
:math:`\gamma^k`, the next iterate of type-II AA can be written as

.. math::
  \begin{align}
  x^{k+1}&=f(x^k)-\sum_{i=0}^{m_k-1}\gamma_i^k\left(f(x^{k-m_k+i+1})- f(x^{k-m_k+i})\right)\\
  &=x^k-g_k-(S_k-Y_k)\gamma^k\\
  &=x^k-(I+(S_k-Y_k)(Y_k^\top Y_k)^{-1}Y_k^\top )g_k\\
  &=x^k-B_kg_k,
  \end{align}

where :math:`S_k=[s_{k-m_k}~\dots~s_{k-1}]`, :math:`s_i=x^{i+1}-x^i` for each
:math:`i`, and

.. math::
   B_k=I+(S_k-Y_k)(Y_k^\top Y_k)^{-1}Y_k^\top

Observe that :math:`B_k` minimizes :math:`\|B_k-I\|_F` subject to
the inverse multi-secant condition :math:`B_kY_k=S_k`, and hence can be regarded
as an approximate inverse Jacobian of :math:`g`. The update of :math:`x^k` can
then be considered as a quasi-Newton-type update, with :math:`B_k` being
a generalized second (or type-II) Broyden's update of :math:`I` satisfying
the inverse multi-secant condition.

Type-I AA
"""""""""

In the same spirit, we define type-I AA, in which we find an approximate
Jacobian of :math:`g` minimizing :math:`\|H_k-I\|_F` subject to the multi-secant
condition :math:`H_kS_k=Y_k`. Assuming that :math:`S_k` is full column rank, we
obtain (by symmetry) that

.. math::
  H_k=I+(Y_k-S_k)(S_k^\top S_k)^{-1}S_k^\top

and the update scheme is defined as

.. math::
  x^{k+1}=x^k-H_k^{-1}g_k

assuming :math:`H_k` to be invertible. A direct application of the Woodbury
matrix identity shows that

.. math::
  H_k^{-1}=I+(S_k-Y_k)(S_k^\top Y_k)^{-1}S_k^\top

where again we have assumed that :math:`S_k^\top Y_k` is invertible.  Notice
that this explicit formula of :math:`H_k^{-1}` is preferred in that the most
costly step, inversion, is implemented only on a small :math:`m_k\times m_k`
matrix.

In SCS
------

In SCS both types of acceleration are available, though by default type-I is
used since it tends to have better performance.  If you wish to use AA then set
the :code:`acceleration_lookback` setting to a non-zero value (10 works well for
many problems and is the default). This setting corresponds to :math:`m`, the
maximum number of SCS iterates that AA will use to extrapolate to the new point.

To enable type-II acceleration then set :code:`acceleration_lookback` to a
negative value, the sign is interpreted as switching the AA type (this is mostly
so that we can test it without fully exposing it the user).

The setting :code:`acceleration_interval` controls how frequently AA is applied.
If :code:`acceleration_interval` :math:`=k` for some integer :math:`k \geq 1`
then AA is applied every :math:`k` iterations (AA simply ignores the
intermediate iterations). This has the benefit of making AA :math:`k` times
faster and approximating a :math:`k` times larger memory, as well as improving
numerical stability by 'decorrelating' the data. On the other hand, older
iterates might be stale.  More work is needed to determine the optimal setting
for this parameter, but 10 appears to work well in practice and is the default.

The details about how the linear systems are solved and updated is abstracted
away into the AA package (eg, QR decomposition, SVD decomposition etc). Exactly
how best to solve and update the equations is still open.

Regularization
""""""""""""""

By default we also add a small amount of regularization to the matrices
that are being inverted in the above expressions, ie, in the type-II update

.. math::
   (Y_k^\top Y_k)^{-1} \text{ becomes } (Y_k^\top Y_k + \epsilon I)^{-1}

for some small :math:`\epsilon > 0`, and similarly for the type-I update

.. math::
   (S_k^\top Y_k)^{-1} \text{ becomes } (S_k^\top Y_k + \epsilon I)^{-1}

which is equivalent to adding regularization to the :math:`S_k^\top S_k` matrix
before using the Woodbury matrix identity.  The regularization ensures the
matrices are invertible and helps stability. In practice type-I tends to require
more regularization than type-II for good performance. The regularization
shrinks the AA update towards the update without AA, since if
:math:`\epsilon\rightarrow\infty` then :math:`\gamma^\star = 0` and the AA step
reduces to :math:`x^{k+1} = f(x^k)`. Note that the regularization can be folded
into the matrices by appending :math:`\sqrt{\epsilon} I` to the bottom of
:math:`S_k` or :math:`Y_k`, which is useful when using a QR or SVD decomposition
to solve the equations.

Max :math:`\gamma` norm
"""""""""""""""""""""""
As the algorithm converges to the fixed point the matrices to be inverted
can become ill-conditioned and AA can become unstable. In this case the
:math:`\gamma` vector can become very large. As a simple heuristic we reject
the AA update and reset the AA state whenever :math:`\|\gamma\|_2` is greater
than :code:`max_weight_norm` (eg, something very large like :math:`10^{10}`).


Safeguarding
""""""""""""

We also apply a safeguarding step to the output of the AA step. Explicitly, let
:math:`x^k` be the current iteration and let :math:`x_\mathrm{AA} = x^{k+1}` be
the output of AA. We reject the AA step if

.. math::
  \|x_\mathrm{AA} - f(x_\mathrm{AA}) \|_2 > \zeta \|x^k - f(x^k) \|_2

where :math:`\zeta` is the safeguarding tolerance factor
(:code:`safeguard_factor`) and defaults to 1. In other words we reject the step
if the norm of the residual after the AA step is larger than some amount (eg, if
it increases the residual from the previous iterate).  After rejecting a step we
revert the iterate to :math:`x^k` and reset the AA state.

Relaxation
""""""""""

In some works relaxation has been shown to improve performance. Relaxation
replaces the final step of AA by mixing the map inputs and outputs as follows:

.. math::
  x^{k+1} = \beta \sum_{j=0}^{m_k}\alpha_j^k f(x^{k-m_k+j}) + (1-\beta) \sum_{j=0}^{m_k}\alpha_j^k x^{k-m_k+j}

where :math:`\beta` is the :code:`relaxation` parameter, and :math:`\beta=1`
recovers vanilla AA. This can be computed using the matrices defined above using

.. math::
  x^{k+1} = \beta (f(x^k) - (S_k - Y_k) \gamma^k) + (1-\beta) (x^k - S_k \gamma^k)


Anderson acceleration API
-------------------------

For completeness, we document the full Anderson acceleration API below.

.. doxygenfile:: include/aa.h
