.. _scaling:

Scaling
=======

In this note we derive the update equations when using a non-identity
scaling when solving :math:`LCP(M, q)`. Standard Douglas-Rachford splitting is:

.. math::
  \begin{align}
  \tilde u^{k+1} &= (I + F)^{-1} w^k \\
  u^{k+1} &= (I + G)^{-1} (2 \tilde u^{k+1} - w^k) \\
  w^{k+1} &= w^k + \alpha (u^{k+1} - \tilde u^{k+1}) \\
  \end{align}


yielding :math:`w^k \rightarrow u^\star + F(u^\star)` where :math:`0 \in
F(u^\star) + G(u^\star)` (assuming such a :math:`u^\star` exists).  Now consider
modifying DR splitting to use a 
diagonal matrix :math:`R` instead of :math:`I`. This is useful because the
matrix :math:`R` can be selected to provide better convergence in practice.
The above becomes:

.. math::
  \begin{align}
  \tilde u^{k+1} &= (R + F)^{-1} R w^k \\
  u^{k+1} &= (R + G)^{-1} R (2 \tilde u^{k+1} - w^k) \\
  w^{k+1} &= w^k + \alpha (u^{k+1} - \tilde u^{k+1}) \\
  \end{align}

which yields :math:`w^k \rightarrow u^\star + R^{-1} F(u^\star)` where :math:`0
\in F(u^\star) + G(u^\star)`.

Cone projection
---------------

In SCS, :math:`F = \mathcal{Q}` and :math:`G = I_{\mathcal{C}_+}`. Note that
:math:`R` has to be chosen so that the cone projection is preserved. To do this
we ensure that the entries of :math:`R > 0` corresponding to each cone are
constant within that cone. This means that R does not affect the cone projection
in any way, since for sub-cone :math:`\mathcal{K}`:

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
described in :ref:`Updating the scale`. Finally,  :math:`d` is determined by
:code:`TAU_FACTOR` in the code defined in :code:`glbopts.h`.



Linear system solve
-------------------

At each iteration SCS must solve the following linear equation:

.. math::
  z = p^k - r \tau, 

where

.. math::
  \begin{align}
  p^k &= (R + M)^{-1} R \mu^k \\
  r   &= (R + M)^{-1} q
  \end{align}

(:math:`R` does *not* appear before :math:`q` in the second expression above).
Now consider :math:`r = (R + M)^{-1} q` and recall 

.. math::
  M = \begin{bmatrix} 
        P  &  A^\top \\
        -A &  0   \\
      \end{bmatrix}


Denote by :math:`R_x = \rho_x I_n` and :math:`R_y = \mathrm{diag}(\rho_y)`.
We want to solve 

.. math::

  \begin{bmatrix}
  r_x \\
  r_y
  \end{bmatrix}
  =
  \begin{bmatrix} 
  R_x + P  &  A^\top \\
  -A &  R_y   \\
  \end{bmatrix}
  \begin{bmatrix}
  q_x \\
  q_y
  \end{bmatrix}

which is quasidefinite if we negate the bottom row:

.. math::

  \begin{bmatrix}
  r_x \\
  r_y
  \end{bmatrix}
  =
  \begin{bmatrix} 
  R_x + P  &  A^\top \\
  A &  -R_y   \\
  \end{bmatrix}
  \begin{bmatrix}
  q_x \\
  -q_y
  \end{bmatrix}

A direct method factorizes the above matrix.
An indirect method can solve via:

.. math::

  \begin{align}
  (R_x + P + A^\top R_y^{-1} A) r_x & = q_x - A^\top R_y^{-1} q_y \\
                            r_y & = R_y^{-1}(A z_x + q_y).
  \end{align}

Vector :math:`v^k`
------------------

In order to get the :math:`v^k` vector that contains :math:`s` and
:math:`\kappa` we use:

.. math::
  v^{k+1} = R( u^{k+1} + w^k - 2 \tilde u^{k+1} ) \rightarrow \mathcal{Q}(u^\star),

and we have

.. math::

  \begin{align}
  v^{k+1} &= R( u^{k+1} + w^k - 2 \tilde u^{k+1} ) \\
          &= R( \Pi_{\mathcal{C}_+} (2 \tilde u^{k+1} - w^k) + w^k - 2 \tilde u^{k+1}) \\
          &= R( \Pi_{\mathcal{C}^*_+} (-2 \tilde u^{k+1} + w^k)) \\
          &= R \tilde v^{k+1} \\
          &\in \mathcal{C}^*_+
  \end{align}

by Moreau and the fact that :math:`R` is chosen to be constant
within each sub-cone :math:`\mathcal{K}`. Finally note that :math:`v^k \perp u^k`, since

.. math::
  \begin{align}
  (u^k)^\top v^k &= (u^k)^\top R \tilde v^{k+1}  \\
   &= \sum_{\mathcal{K} \in \mathcal{C}_+} r_\mathcal{K} (u^k_\mathcal{K})^\top \tilde v^k_\mathcal{K} \\
   &= 0
  \end{align}

since :math:`\tilde v^k \perp u^k` and :math:`R` is chosen to be constant
within each sub-cone :math:`\mathcal{K}`.

Root-plus function
------------------

Finally, the :code:`root_plus` function is modified to be the solution
of the following quadratic equation:

.. math::
  \tau^2 (d + r^\top R r) + \tau (r^\top R \mu^k - 2 r^\top R p^k - d \eta^k) + p^k R (p^k - \mu^k) = 0.

Other than when computing :math:`\kappa` (which does not affect the algorithm)
this is the *only* place where :math:`d` appears, so we have a lot of
flexibility in how to choose it and it can even change from iteration to
iteration. It is an open question on how best to select this parameter.  See the
:code:`dot_with_diag_scaling` function in :code:`src/scs.c`.

