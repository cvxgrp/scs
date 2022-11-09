.. _scaling:

Non-identity DR Scaling
=======================

In this note we derive the update equations when using a non-identity scaling.
Standard Douglas-Rachford splitting applied to the SCS problem is described in
the :ref:`algorithm` page.  Now consider modifying DR splitting to use a
diagonal matrix :math:`R` instead of :math:`I`. This is useful because the
matrix :math:`R` can be selected to provide better convergence in practice.  The
algorithm then becomes:

.. math::
  \begin{align}
  \tilde u^{k+1} &= (R + \mathcal{Q})^{-1} R w^k \\
  u^{k+1} &= (R + N_{\mathcal{C}_+})^{-1} R (2 \tilde u^{k+1} - w^k) \\
  w^{k+1} &= w^k + u^{k+1} - \tilde u^{k+1} \\
  \end{align}

which yields :math:`w^k \rightarrow u^\star + R^{-1} \mathcal{Q}(u^\star)` where
:math:`0 \in \mathcal{Q}(u^\star) + N_{\mathcal{C}_+}(u^\star)`.

This changes the first two steps of the procedure. The :ref:`linear projection
<linear_solver>` and the cone projection (explained next).

Cone projection
---------------

Note that :math:`R` has to be chosen so that the cone projection is preserved.
To do this we ensure that the entries of :math:`R > 0` corresponding to each
cone are constant within that cone (with the box cone being a slight exception
to this). This means that R does not affect the cone projection in any way,
since for sub-cone :math:`\mathcal{K}`:

.. math::
   \begin{align}
    &x = (r I + \partial I_{\mathcal{K}})^{-1} (r y) \\
    \Rightarrow \quad  &  0 \in r(x - y) + \partial I_{\mathcal{K}}(x) \\
    \Rightarrow \quad & x = \mbox{argmin}_z || z - y ||^2 \mbox{ s.t. } z \in \mathcal{K}.
    \end{align}

In other words, :math:`R` is selected such that

.. math::
   (R +  \partial I_{\mathcal{C}_+})^{-1} R = (I +  \partial I_{\mathcal{C}_+})^{-1}


Selecting :math:`R`
-------------------

For SCS we take

.. math::
  R =  \begin{bmatrix}
    \rho_x I_n   &        0    &   0 \\
    0     &     \mathrm{diag}(\rho_y) &   0 \\
    0     &               0     &   d
  \end{bmatrix}

where :math:`I_n` is :math:`n \times n` identity, :math:`\rho_x \in \mathbf{R}`,
:math:`\rho_y \in \mathbf{R}^m` and :math:`d \in \mathbf{R}`. The :math:`\rho_y`
term includes the effect of the parameter :code:`scale`, which is updated
heuristically to improve convergence. Basically

.. math::

    \rho_y \approx (1/\mathrm{scale}) I

with the only difference being that each cone can modify this
relationship slightly (but currently only the zero cone does, which sets
:math:`\rho_y = 1 / (1000\ \mathrm{scale})`).
So as scale increases, :math:`\rho_y` decreases and vice-versa.

The quantity :math:`\rho_x` is determined by the :ref:`setting <settings>` value
:code:`rho_x`. The quantity :math:`\rho_y` is determined by the :ref:`setting
<settings>` value :code:`scale` but is updated adaptively using the techniques
described in :ref:`Dynamic scale updating <updating_scale>`. Finally,  :math:`d`
is determined by :code:`TAU_FACTOR` in the code defined in :code:`glbopts.h`.

Root-plus function
------------------

Finally, the :code:`root_plus` function is modified to be the solution
of the following quadratic equation:

.. math::
  \tau^2 (d + r^\top R_{-1} r) + \tau (r^\top R_{-1} \mu^k - 2 r^\top R_{-1} p^k - d \eta^k) + p^k R_{-1} (p^k - \mu^k) = 0,

where :math:`R_{-1}` corresponds to the first :math:`n+m` entries of :math:`R`.
Other than when computing :math:`\kappa` (which does not affect the algorithm)
this is the *only* place where :math:`d` appears, so we have a lot of
flexibility in how to choose it and it can even change from iteration to
iteration. It is an open question on how best to select this parameter.

Moreau decomposition
--------------------
The projection onto a cone under a diagonal scaling also satisfies a
Moreau-style decomposition identity, as follows:

.. math::
   x + R^{-1} \Pi_\mathcal{C^*}^{R^{-1}} ( - R x ) = \Pi_\mathcal{C}^R ( x )

where :math:`\Pi_\mathcal{C}^R` denotes the projection onto convex cone
:math:`\mathcal{C}` under the :math:`R`-norm, which is defined as

.. math::
  \|x\|_R = \sqrt{x^\top R x}.

This identity is useful when deriving cone projection routines, though most
cones are either invariant to this or we enforce that :math:`R` is constant
across them. Note that the two components of the decomposition are
:math:`R`-orthogonal.

Dual vector
-----------

In order to get the dual vector :math:`v^k` that contains :math:`s^k` and
:math:`\kappa^k` we use:

.. math::
  v^{k+1} = R( u^{k+1} + w^k - 2 \tilde u^{k+1} ) \rightarrow \mathcal{Q}(u^\star),

and we have

.. math::
  \begin{align}
  v^{k+1} &= R( u^{k+1} + w^k - 2 \tilde u^{k+1} ) \\
          &= R( \Pi^R_{\mathcal{C}_+} (2 \tilde u^{k+1} - w^k) + w^k - 2 \tilde u^{k+1}) \\
          &= R( R^{-1} \Pi^{R^{-1}}_{\mathcal{C}^*_+} (R(w^k -2 \tilde u^{k+1}))) \\
          &= \Pi^{R^{-1}}_{\mathcal{C}^*_+} (R(w^k -2 \tilde u^{k+1})) \\
          &\in \mathcal{C}^*_+
  \end{align}

by Moreau, and finally note that :math:`v^k \perp
u^k` from the fact that the Moreau decomposition is :math:`R`-orthogonal.

.. _updating_scale:

Dynamic scale updating
----------------------
The choice of the :code:`scale` parameter can have a large impact on the
performance of the algorithm and the optimal choice is highly problem
dependent. SCS can dynamically adjust the :code:`scale` parameter
on the fly via a heuristic procedure that can substantially improve convergence
in practice. This procedure is enabled by the :code:`adaptive_scale`
:ref:`setting <settings>`. The procedure attempts to balance the convergence
rate of the primal residual with the dual residual. Loosely speaking, the
:code:`scale` parameter will be increased if the primal residual is much larger
than the dual and decreased if the opposite is true.

Specifically, at iteration :math:`k` consider the case where :math:`l`
iterations have elapsed since the last update of the :code:`scale` parameter,
and denote by :math:`(x, y, \tau) = u^k` and :math:`(0, s, \kappa) = v^k`, and
the *relative* residuals as

.. math::
   \hat r^k_p = \frac{\|Ax + s - b \tau\|}{\max(\|Ax\|, \|s\|, \|b \tau \|)}

.. math::
   \hat r^k_d = \frac{\|Px + A^\top y + c \tau\|}{\max(\|Px\|, \|A^\top y\|, \|c \tau \|)}

where by default we use the :math:`\ell_\infty` norm for these quantities,
but can be changed using the :code:`SCALE_NORM` constant in
:code:`include/glbopts.h`.
Now consider

.. math::
  \beta = \left(\prod_{i=0}^{l-1} \frac{\hat r^{k-i}_p}{\hat r^{k-i}_d}\right)^{1/l}

ie, :math:`\beta` corresponds to the geometric mean of the ratio of the relative
residuals across the last :math:`l` iterations. If this number is larger than a
constant (eg, 3) or smaller than another constant (eg, 1/3) *and* if sufficient
iterations have passed since the last update (eg, 100, as determined by
:code:`RESCALING_MIN_ITERS`) then an update of the :code:`scale` parameter is
triggered:

.. math::
   \mbox{scale}^+ = \sqrt{\beta}\ \mbox{scale}

The presence of the square root is to prevent over-shooting the 'optimal'
scale parameter, which could lead to oscillation.

Note that if the :ref:`linear system <linear_solver>` is being solved using a
direct method, then updating the :code:`scale` parameter will require a new
factorization of the perturbed matrix, so is somewhat expensive for larger
problems and should be done sparingly. Also, since the changing the
:code:`scale` changes the operator we are using in DR splitting we also need to
perform a reset of the :ref:`Anderson acceleration <acceleration>`.

