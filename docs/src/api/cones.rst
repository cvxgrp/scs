.. _cones:

Cones
-----

The cone :math:`\mathcal{K}` can be any Cartesian product of the following primitive cones.


.. list-table::
   :widths: 30 30 40
   :header-rows: 1

   * - Name
     - Description
     - Entries in :ref:`ScsCone <ScsCone>`

   * - Zero cone
     - :math:`\{s \mid s = 0 \}`
     - :code:`z` length of cone.
   * - Positive orthant
     - :math:`\{s \mid s \geq 0 \}`
     - :code:`l` length of cone.
   * - Box cone
     - :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^n \mid t l \leq s \leq t u  \}`
     - :code:`bu, bl` correspond to :math:`u,l` and are arrays of len :code:`bsize` :math:`=n`. If :code:`bsize=0` cone is disabled.
   * - Second-order cone
     - :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^n\mid \|s\|_2 \leq t  \}`
     - :code:`q` array of SOC lengths with :code:`qsize` elements, each :math:`q[i] = n_i+1`.
   * - Positive semidefinite cone
     - :math:`\{ S \in \mathbf{R}^{n \times n} \mid \min_i \lambda_i(S) \geq 0, S = S^\top  \}`
     - :code:`s` array of PSD cone lengths with :code:`ssize` elements, each :math:`s[i] = n_i`. See :ref:`note <sdcone>` about scaling.
   * - Exponential cone
     - :math:`\{   (x,y,z) \in \mathbf{R}^3 \mid y e^{x/y} \leq z, y>0  \}`
     - :code:`ep` number of cone triples.
   * - Dual exponential cone
     - :math:`\{  (u,v,w)\in \mathbf{R}^3 \mid -u e^{v/u} \leq e w, u<0 \}`
     - :code:`ed` number of cone triples.
   * - Power cone
     - :math:`\{  (x,y,z) \in \mathbf{R}^3 \mid x^p y^{1-p} \geq |z|\}`
     - :code:`p` array of :math:`p\in[-1,1]` powers with :code:`psize` elements, positive entries correspond to power cone.
   * - Dual power cone
     - :math:`\{ (u,v,w)\in \mathbf{R}^3 \mid \left(\frac{u}{p}\right)^p \left(\frac{v}{1-p}\right)^{1-p} \geq |w|\}`
     - :code:`p` array :math:`p\in[-1,1]` powers with :code:`psize` elements, negative entries correspond to dual power cone (sign is flipped).


**Note**:
The rows of the data matrix :math:`A` correspond to the cones in :math:`\mathcal{K}`. The rows of
:math:`A` must be in the order of the cones given above, i.e., first come the rows that
correspond to the zero cones, then those that correspond to the positive
orthants, then the box cone, etc.

.. _sdcone:

Semidefinite cone
^^^^^^^^^^^^^^^^^^^^^^^^

SCS assumes that the matrix variables and the input data corresponding to
semidefinite cones have been vectorized by scaling the off-diagonal entries by
:math:`\sqrt{2}` and stacking the lower triangular elements column-wise. For a :math:`k \times k`
matrix variable (or data matrix) this operation would create a vector of length
:math:`k(k+1)/2`. Scaling by :math:`\sqrt{2}` is required to preserve the inner-product.

To recover the matrix solution this operation must be inverted on the components
of the vector returned by SCS corresponding to semidefinite cones. That is, the
off-diagonal entries must be scaled by :math:`1/\sqrt{2}` and the upper triangular
entries are filled in by copying the values of lower triangular entries.

More explicitly, we want to express :math:`Tr(C X)` as :math:`\text{vec}(C)^\top \text{vec}(X)`,
where the :math:`\text{vec}` operation takes the :math:`k \times k` matrix

.. math::

  X =  \begin{bmatrix}
          X_{11} & X_{12} & \ldots & X_{1k}  \\
          X_{21} & X_{22} & \ldots & X_{2k}  \\
          \vdots & \vdots & \ddots & \vdots  \\
          X_{k1} & X_{k2} & \ldots & X_{kk}  \\
        \end{bmatrix}

and produces a vector consisting of the lower triangular elements scaled and arranged as

.. math::

  \text{vec}(X) = (X_{11}, \sqrt{2} X_{21}, ..., \sqrt{2} X_{k1}, X_{22}, \sqrt{2}X_{32}, ..., \sqrt{2}X_{k(k-1)}, X_{kk}).

