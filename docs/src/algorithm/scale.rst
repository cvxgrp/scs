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

Note that :math:`R` has to be chosen so that the cone projection can still be
computed. To do this we ensure that the entries of :math:`R > 0` corresponding
to each cone are constant within that cone, with two exceptions. Both arise
only with the :code:`adaptive_diag_scale` :ref:`setting <settings>`, which is
on by default and is described in :ref:`per_row_refinement`: each row of the
box cone can have its own entry, and the entries of a real PSD block can take
a rank-one form. (The zero and nonnegative cones are products of
one-dimensional cones, so each of their rows can have its own entry too.)
When the entries are constant within a cone, :math:`R` does not affect the
projection onto that cone in any way, since for sub-cone :math:`\mathcal{K}`:

.. math::
    &x = (r I + \partial I_{\mathcal{K}})^{-1} (r y) \\
    \Rightarrow \quad  &  0 \in r(x - y) + \partial I_{\mathcal{K}}(x) \\
    \Rightarrow \quad & x = \mbox{argmin}_z || z - y ||^2 \mbox{ s.t. } z \in \mathcal{K}.

In other words, for these cones

.. math::
   (R +  \partial I_{\mathcal{K}})^{-1} R = (I +  \partial I_{\mathcal{K}})^{-1}.

For the two exceptions the projection does depend on :math:`R`, and it is
computed exactly in the norm that :math:`R` defines (see the Moreau
decomposition below). The box cone projection solves a weighted problem in
its scalar variable, and the PSD projection applies a diagonal congruence
before and after the usual eigendecomposition.


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

with two differences. First, each cone can modify this relationship slightly
(currently only the zero cone does, which sets
:math:`\rho_y = 1 / (1000\ \mathrm{scale})`). Second, with
:code:`adaptive_diag_scale` each row :math:`i` has its own multiplier
:math:`f_i`, which acts like a scale for that row alone:

.. math::

    \rho_{y,i} = \frac{1}{\mathrm{scale}\ f_i}

for a row outside the zero cone, and
:math:`\rho_{y,i} = 1 / (1000\ \mathrm{scale}\ f_i)` for a row of the zero
cone. The multipliers start at one, so at the start of a solve on a new
workspace they change nothing, and they are kept in
:math:`[10^{-3}, 10^{3}]`. How they are updated is described in
:ref:`per_row_refinement`. So as scale increases, :math:`\rho_y` decreases
and vice-versa.

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

The procedure does not run at every iteration. It runs at the iterations where
SCS computes the residuals to test for convergence, which is every 25
iterations starting at iteration 0 (:code:`CONVERGED_INTERVAL` in
:code:`include/glbopts.h`). At each of these iterations it takes one sample
of the imbalance between the two residuals and then decides whether to update.

Specifically, at such an iteration :math:`k` consider the case where :math:`l`
iterations have elapsed since the last update of the :code:`scale` parameter
or of the per-row multipliers described below, and denote by
:math:`(x, y, \tau) = u^k` and :math:`(0, s, \kappa) = v^k`.
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
   \theta^k = (1 - \sigma^k_p)(1 - \sigma^k_d),

so that the scalar :math:`\theta^k` is close to one when both residuals are
small against their terms and zero when either is at its ceiling. The sample
taken at iteration :math:`k` is

.. math::
   \rho^k = \left(\frac{d^k_p}{d^k_d}\right)^{\theta^k}
            \left(\frac{\|\tilde r^k_p\|}{\|\tilde r^k_d\|}\right)^{1 - \theta^k}.

Since the stopping thresholds enter the sample, the :code:`eps_abs` and
:code:`eps_rel` settings affect the updates of the :code:`scale`, and
therefore the iterates, not only the point at which the solve stops. They
enter only through the ratio
:math:`\epsilon_\mathrm{abs} / \epsilon_\mathrm{rel}`: a common factor cancels
in :math:`d^k_p / d^k_d` and in the saturations. So tightening both
tolerances by the same factor, for example from :math:`10^{-4}` to
:math:`10^{-6}` each, leaves the samples unchanged up to rounding, while
changing one of them alone changes the samples, and the solve can then follow
a different path from the first update on.

Now let :math:`S` be the set of iterations at which a sample has been taken
since the last update (or since the start of the solve), up to and including
:math:`k`, and consider

.. math::
  \beta = \left(\prod_{j \in S} \rho^{j}\right)^{1/|S|}

ie, :math:`\beta` is the geometric mean of the samples taken since the last
update. If this number is larger than 10 or smaller than 1/10 *and* if
sufficient iterations have passed since the last update (:math:`l` at least
100, as determined by :code:`RESCALING_MIN_ITERS`) then an update of the
:code:`scale` parameter is triggered:

.. math::
   \mbox{scale}^+ = \min(\max(\sqrt{\beta}\ \mbox{scale},\ 10^{-6}),\ 10^{6}).

The presence of the square root is to prevent over-shooting the 'optimal'
scale parameter, which could lead to oscillation. The limits are
:code:`MIN_SCALE_VALUE` and :code:`MAX_SCALE_VALUE` in
:code:`include/glbopts.h`. If the :code:`scale` is already at one of them and
:math:`\beta` would push it further out, then no update of the :code:`scale`
takes place (the per-row multipliers described below can still be updated on
their own). The limits apply to the updates only: the :code:`scale`
:ref:`setting <settings>` can be any positive number.

With these values :math:`\beta` is the geometric mean of four samples the
first time an update is allowed after an earlier one (five for the first
update of a solve, which also has the sample from iteration 0), and it gains
one sample every 25 iterations for as long as no update is triggered.

Note that if the :ref:`linear system <linear_solver>` is being solved using a
direct method, then every update requires a new factorization of the
perturbed matrix, so is somewhat expensive for larger problems and should be
done sparingly. Also, since changing the :code:`scale` or the per-row
multipliers changes the operator we are using in DR splitting, every update
also resets the :ref:`Anderson acceleration <acceleration>`, whichever linear
solver is in use. After a reset the acceleration has to fill its memory again
before it produces a step: the first call after the reset only records the
iterate, and the first step comes :code:`acceleration_lookback` calls later,
with one call every :code:`acceleration_interval` iterations. With the default
settings that is 55 iterations without an accelerated step after each update,
and updates can be as close as 100 iterations apart.

The number of updates made in a solve is returned in :code:`scale_updates` in
the :ref:`info <info>` struct, and the final :code:`scale` in :code:`scale`.


.. _per_row_refinement:

Per-row refinement
------------------

The :code:`scale` parameter above moves the whole :math:`R_y` block together.
With :code:`adaptive_diag_scale` each row :math:`i` carries its own multiplier
:math:`f_i` as well. It starts at one and divides that row's entry of
:math:`\rho_y`, so row :math:`i` behaves as if its scale were
:math:`f_i\ \mathrm{scale}`.

Both :code:`adaptive_scale` and :code:`adaptive_diag_scale` are on by default.
The multipliers are computed inside the procedure of the previous section,
which only runs when :code:`adaptive_scale` is on, so
:code:`adaptive_diag_scale` requires :code:`adaptive_scale`. If
:code:`adaptive_scale` is off then :code:`adaptive_diag_scale` is switched off
as well, without a warning, and :math:`R` stays fixed for the whole solve. If
:code:`adaptive_scale` is on and :code:`adaptive_diag_scale` is off then only
the scalar :code:`scale` is updated.

The multipliers are updated from each row's relative primal residual

.. math::
   \hat r^k_{p,i} = \frac{|(Ax + s - b\tau)_i|}
                         {\max(|(Ax)_i|, |s_i|, |b_i \tau|, \epsilon_{\rm den})},

by the damped multiplicative step

.. math::
   f_i^+ = f_i \left(\hat r^k_{p,i} / g\right)^{1/4},

where :math:`g` is the geometric mean of :math:`\hat r^k_{p,i}` over all
:math:`m` rows. A row whose relative residual is above the geometric mean
therefore gets a larger multiplier (a larger scale for that row), and a row
below it a smaller one. The result is clipped to :math:`[10^{-3}, 10^{3}]`
(the exponent and the range are :code:`DIAG_SCALE_DAMP`,
:code:`DIAG_SCALE_MULT_MIN` and :code:`DIAG_SCALE_MULT_MAX` in
:code:`include/glbopts.h`) and then constrained within each non-polyhedral
cone block to a form the cone projection can absorb. For most cone types that
is a single value, the arithmetic mean of the block's clipped multipliers. For real PSD
blocks it is the richer rank-one family :math:`w_{ij} = \delta_i \delta_j`
(a diagonal congruence :math:`X \to D X D`, which maps the cone to itself),
fitted per block by least squares in log space, with each :math:`\delta_i`
clipped so that every :math:`w_{ij}` stays in the same range. The PSD
projection under such a metric is computed exactly via a congruence sandwich
around the eigendecomposition. This is the dynamic analog of the row half of
Ruiz equilibration, using the residuals the algorithm actually produces rather
than the matrix norms.

The multipliers are updated whenever the scalar :code:`scale` is. They are
also updated on their own, with a new factorization and a reset of the
acceleration exactly as for the scalar, when some multiplier would move by
more than the threshold used for the scalar (a factor of :math:`\sqrt{10}`)
and at least 100 iterations have passed since the last update of either
kind, so that a scalar :code:`scale` sitting at its limit cannot freeze the
diagonal. The move is measured on the multipliers as they would be
installed, after each cone block has been given its form. A single row of a
second-order or PSD block that asks for a large step therefore does not cost a
factorization and a reset when the value its block ends up with has barely
changed. Rows of the zero, nonnegative and box cones each keep their own
multiplier, so for them nothing is averaged away. An update of the multipliers
on their own is counted in :code:`scale_updates` (see :ref:`info`), but the
:code:`scale` printed in the log and returned in the info struct does not
change.

The floor :math:`\epsilon_{\rm den}` in the denominator is :math:`10^{-3}`
(:code:`DEN_FLOOR_FRAC` in :code:`include/glbopts.h`) times the
root-mean-square over all rows of the denominator without the floor,
:math:`\max(|(Ax)_i|, |s_i|, |b_i \tau|)`. It is recomputed from the current
iterate each time, and the rows are taken in the space the iteration runs in
(that is, after :ref:`equilibration <equilibration>` when :code:`normalize` is
on), like :math:`\tilde r^k_p` above. It matters because a row
whose :math:`b_i` is zero divides only by :math:`\max(|(Ax)_i|, |s_i|)`, and
both of those shrink as the iterates converge, so that row's ratio inflates
even when it is converging perfectly well. Without the floor the profile
partly reports which terms a row happens to contain rather than how the row is
converging, and the metric chases that artifact. Rows whose denominators are
of ordinary size are unaffected.
