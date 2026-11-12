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
  \tilde u^{k+1} &= (R + \mathcal{Q})^{-1} R w^k \\
  u^{k+1} &= (R + N_{\mathcal{C}_+})^{-1} R (2 \tilde u^{k+1} - w^k) \\
  w^{k+1} &= w^k + u^{k+1} - \tilde u^{k+1} \\

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
    &x = (r I + \partial I_{\mathcal{K}})^{-1} (r y) \\
    \Rightarrow \quad  &  0 \in r(x - y) + \partial I_{\mathcal{K}}(x) \\
    \Rightarrow \quad & x = \mbox{argmin}_z || z - y ||^2 \mbox{ s.t. } z \in \mathcal{K}.

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
  v^{k+1} &= R( u^{k+1} + w^k - 2 \tilde u^{k+1} ) \\
          &= R( \Pi^R_{\mathcal{C}_+} (2 \tilde u^{k+1} - w^k) + w^k - 2 \tilde u^{k+1}) \\
          &= R( R^{-1} \Pi^{R^{-1}}_{\mathcal{C}^*_+} (R(w^k -2 \tilde u^{k+1}))) \\
          &= \Pi^{R^{-1}}_{\mathcal{C}^*_+} (R(w^k -2 \tilde u^{k+1})) \\
          &\in \mathcal{C}^*_+

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
and denote by :math:`(x, y, \tau) = u^k` and :math:`(0, s, \kappa) = v^k`.
Two measures of the imbalance between the primal and the dual residual are
blended. The first is how far each residual is from its own :ref:`stopping
threshold <termination>`, with :math:`\tau` multiplied through:

.. math::
   d^k_p = \frac{\|Ax + s - b \tau\|}{\epsilon_\mathrm{abs} \tau + \epsilon_\mathrm{rel} \max(\|Ax\|, \|s\|, \|b \tau \|)}

.. math::
   d^k_d = \frac{\|Px + A^\top y + c \tau\|}{\epsilon_\mathrm{abs} \tau + \epsilon_\mathrm{rel} \max(\|Px\|, \|A^\top y\|, \|c \tau \|)}

so a value of one or less means that residual would pass the test. By default
we use the :math:`\ell_\infty` norm for these quantities, which can be changed
using the :code:`SCALE_NORM` constant in :code:`include/glbopts.h`. The second
measure is the residuals themselves, :math:`\tilde r^k_p` and
:math:`\tilde r^k_d`, in the space the iteration runs in (that is, after
:ref:`equilibration <equilibration>` when :code:`normalize` is on), with
nothing divided out.

Balancing the distances :math:`d^k_p` and :math:`d^k_d` evens out how far each
side still has to go, which is the right target whenever both can be trusted.
A distance cannot be trusted when the residual is as large as the terms it is
measured against: such a residual has no correct digits, and its distance then
sits at a ceiling of about :math:`1 / \epsilon_\mathrm{rel}` whatever the
iterate does. The clearest case is a feasibility problem, :math:`c = 0` and
:math:`P = 0`, where the dual residual is :math:`A^\top y` and is therefore
the only term it is measured against. Its distance is at the ceiling on every
iteration, so balancing the distances alone would read the dual as permanently
behind and drive the :code:`scale` down to its lower limit. To handle this,
define the saturation of each side and the resulting weight

.. math::
   \sigma^k_p = \min(1, \epsilon_\mathrm{rel} d^k_p), \quad
   \sigma^k_d = \min(1, \epsilon_\mathrm{rel} d^k_d), \quad
   w^k = (1 - \sigma^k_p)(1 - \sigma^k_d),

so that :math:`w^k` is close to one when both residuals are small against
their terms and zero when either is at its ceiling. The sample taken at
iteration :math:`k` is

.. math::
   \rho^k = \left(\frac{d^k_p}{d^k_d}\right)^{w^k}
            \left(\frac{\|\tilde r^k_p\|}{\|\tilde r^k_d\|}\right)^{1 - w^k}.

Now consider

.. math::
  \beta = \left(\prod_{i=0}^{l-1} \rho^{k-i}\right)^{1/l}

ie, :math:`\beta` corresponds to the geometric mean of the samples
across the last :math:`l` iterations. If this number is larger than a
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


Per-row refinement
------------------

The :code:`scale` parameter above moves the whole :math:`R_y` block together.
With :code:`adaptive_diag_scale` each row carries its own multiplier as well,
updated at scale-update events from that row's relative primal residual

.. math::
   \hat r^k_{p,i} = \frac{|(Ax + s - b\tau)_i|}
                         {\max(|(Ax)_i|, |s_i|, |b_i \tau|, \epsilon_{\rm den})},

by a damped multiplicative step towards the geometric mean of the profile,
clamped to a bounded range and constrained within each non-polyhedral cone
block to a form the cone projection can absorb: a single metric value for
most cone types, but for real PSD blocks the richer rank-one family
:math:`w_{ij} = \delta_i \delta_j` (a diagonal congruence
:math:`X \to D X D`, which maps the cone to itself), fitted per block by
least squares in log space. The PSD projection under such a metric is
computed exactly via a congruence sandwich around the eigendecomposition. An update
also fires when the accumulated per-row drift alone is large enough, so that
a railed scalar :code:`scale` cannot freeze the diagonal. This is the
dynamic analog of the row half of Ruiz equilibration, using the residuals the
algorithm actually produces rather than the matrix norms.

The floor :math:`\epsilon_{\rm den}` in the denominator is a small fraction
of the root-mean-square denominator over all rows. It matters because a row
whose :math:`b_i` is zero divides only by :math:`\max(|(Ax)_i|, |s_i|)`, and
both of those shrink as the iterates converge, so that row's ratio inflates
even when it is converging perfectly well. Without the floor the profile
partly reports which terms a row happens to contain rather than how the row is
converging, and the metric chases that artifact. Rows whose denominators are
of ordinary size are unaffected.
