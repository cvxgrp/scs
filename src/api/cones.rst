.. _cones:

Cones
=====

The cone :math:`\mathcal{K}` can be any Cartesian product of the following primitive cones.


.. list-table::
   :header-rows: 1

   * - Name
     - Description
     - Entries in :ref:`ScsCone <ScsCone>`
     - Number of rows in :math:`A, b`
   * - Zero cone
     - :math:`\{s \mid s = 0 \}`
     - :code:`z` length of cone.
     - :code:`z`
   * - Positive (or linear) cone
     - :math:`\{s \mid s \geq 0 \}`
     - :code:`l` length of cone.
     - :code:`l`
   * - Box cone
     - :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^k \mid t l \leq s \leq t u  \}`
     - :code:`bu, bl` correspond to :math:`u,l` and are arrays of length :math:`\max(\text{bsize}-1, 0)` where :code:`bsize` is total cone length.
     - :code:`bsize`
   * - Second-order cone
     - :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^k \mid \|s\|_2 \leq t  \}`
     - :code:`q` array of SOC lengths with :code:`qsize` elements, each :math:`q_i` is total length of that cone.
     - :math:`\displaystyle \sum_{i=1}^{\text{qsize}} q_i`
   * - Positive semidefinite cone
     - :math:`\{ s \in \mathbf{R}^{k(k+1)/2} \mid \text{mat}(s) \succeq 0 \}` (See :ref:`note <sdcone>`)
     - :code:`s` array of PSD cone lengths with :code:`ssize` elements, each :math:`s_i = k_i` in description.
     - :math:`\displaystyle \sum_{i=1}^{\text{ssize}} \frac{s_i(s_i+1)}{2}`
   * - Exponential cone
     - :math:`\{   (x,y,z) \in \mathbf{R}^3 \mid y e^{x/y} \leq z, y>0  \}`
     - :code:`ep` number of cone triples.
     - :math:`3 \times`:code:`ep`
   * - Dual exponential cone
     - :math:`\{  (u,v,w)\in \mathbf{R}^3 \mid -u e^{v/u} \leq e w, u<0 \}`
     - :code:`ed` number of cone triples.
     - :math:`3 \times`:code:`ed`
   * - Power cone
     - :math:`\{  (x,y,z) \in \mathbf{R}^3 \mid x^p y^{1-p} \geq |z|\}`
     - :code:`p` array of :math:`p\in[-1,1]` powers with :code:`psize` elements, positive entries correspond to power cone.
     - :math:`3 \times`:code:`psize` (total for primal and dual power cone)
   * - Dual power cone
     - :math:`\{ (u,v,w)\in \mathbf{R}^3 \mid \left(\frac{u}{p}\right)^p \left(\frac{v}{1-p}\right)^{1-p} \geq |w|\}`
     - :code:`p` array :math:`p\in[-1,1]` powers with :code:`psize` elements, negative entries correspond to dual power cone (sign is flipped).
     - :math:`3 \times`:code:`psize` (total for primal and dual power cone)


**Note**:
The rows of the data matrix :math:`A` correspond to the cones in
:math:`\mathcal{K}`. The rows of :math:`A` **must be in the order of the cones
given above**, i.e., first come the rows that correspond to the zero cones, then
those that correspond to the positive orthants, then the box cone, etc.

Within a cone the rows should be ordered to agree with the mathematical
description of the cone as given above. For instance, the box cone is defined
as :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^k \mid t l \leq s \leq t u
\}` and consequently the variable vector is stacked as :code:`[t;s]`, ie,
:code:`t` comes first then :code:`s`. Similarly, the exponential cone is
defined as :math:`\{   (x,y,z) \in \mathbf{R}^3 \mid y e^{x/y} \leq z, y>0  \}`
and therefore the variable vector is stacked as :code:`[x, y, z]`, etc.

.. _sdcone:

Semidefinite cones
------------------

The symmetric positive semidefinite cone of matrices is the set

.. math::
   \{S \in \mathbf{R}^{k \times k} \mid  S = S^\top,  x^\top S x \geq 0 \ \forall x \in \mathbf{R}^k \}

and for short we use :math:`S \succeq 0` to denote membership. SCS
vectorizes this cone in a special way which we detail here.

SCS assumes that the input data corresponding to
semidefinite cones have been vectorized by scaling the off-diagonal entries by
:math:`\sqrt{2}` and stacking the lower triangular elements column-wise. For a :math:`k \times k`
matrix variable (or data matrix) this operation would create a vector of length
:math:`k(k+1)/2`. Scaling by :math:`\sqrt{2}` is required to preserve the inner-product.

**This must be done for the rows of both** :math:`A` **and** :math:`b` **that correspond to semidefinite cones and must be done independently for each semidefinite cone.**

More explicitly, we want to express :math:`\text{Trace}(Y S)` as :math:`\text{vec}(Y)^\top \text{vec}(S)`,
where the :math:`\text{vec}` operation takes the (assumed to be symmetric) :math:`k \times k` matrix

.. math::

  S =  \begin{bmatrix}
          S_{11} & S_{12} & \ldots & S_{1k}  \\
          S_{21} & S_{22} & \ldots & S_{2k}  \\
          \vdots & \vdots & \ddots & \vdots  \\
          S_{k1} & S_{k2} & \ldots & S_{kk}  \\
        \end{bmatrix}

and produces a vector consisting of the lower triangular elements scaled and arranged as

.. math::

  \text{vec}(S) = (S_{11}, \sqrt{2} S_{21}, \ldots, \sqrt{2} S_{k1}, S_{22}, \sqrt{2}S_{32}, \dots, S_{k-1,k-1}, \sqrt{2}S_{k,k-1}, S_{kk}) \in \mathbf{R}^{k(k+1)/2}.

To recover the matrix solution this operation must be inverted on the components
of the vectors returned by SCS corresponding to each semidefinite cone. That is, the
off-diagonal entries must be scaled by :math:`1/\sqrt{2}` and the upper triangular
entries are filled in by copying the values of lower triangular entries.
Explicitly, the inverse operation takes vector :math:`s \in
\mathbf{R}^{k(k+1)/2}` and produces the matrix

.. math::
  \text{mat}(s) =  \begin{bmatrix}
                    s_{1} & s_{2} / \sqrt{2} & \ldots & s_{k} / \sqrt{2}  \\
                    s_{2} / \sqrt{2} & s_{k+1} & \ldots & s_{2k-1} / \sqrt{2}  \\
                    \vdots & \vdots & \ddots & \vdots  \\
                    s_{k} / \sqrt{2} & s_{2k-1} / \sqrt{2} & \ldots & s_{k(k+1) / 2}  \\
                    \end{bmatrix}
  \in \mathbf{R}^{k \times k}.


So the cone definition that SCS uses is

.. math::
   \mathcal{S}_+^k = \{ \text{vec}(S) \mid S \succeq 0\} = \{s \in \mathbf{R}^{k(k+1)/2} \mid \text{mat}(s) \succeq 0 \}.

Example
^^^^^^^

For a concrete example in python see :ref:`py_mat_completion`.
Here we consider the symmetric positive semidefinite cone constraint over
variables :math:`x \in \mathbf{R}^n` and :math:`S \in \mathbf{R}^{k \times k}`

.. math::
    B - \sum_{i=1}^n \mathcal{A}_i x_i = S \succeq 0

where data :math:`B, \mathcal{A}_1, \ldots, \mathcal{A}_n \in \mathbf{R}^{k
\times k}` are symmetric. We can write this in the canonical form over a new
variable :math:`s \in \mathcal{S}_+^k`:

.. math::
  \begin{align}
  s &= \text{vec}(S)\\
    &= \text{vec}(B - \sum_{i=1}^n \mathcal{A}_i x_i) \\
    &= \text{vec}(B) - \sum_{i=1}^n \text{vec}(\mathcal{A}_i) x_i \\
    &= b - Ax
  \end{align}

using the fact that :math:`\text{vec}` is linear, where :math:`b =
\text{vec}(B)` and

.. math::
  A =
  \begin{bmatrix}
   \text{vec}(\mathcal{A}_1) & \text{vec}(\mathcal{A}_2) & \cdots & \text{vec}(\mathcal{A}_n)
  \end{bmatrix}

i.e., the vectors :math:`\text{vec}(\mathcal{A}_i)` stacked columnwise.
This is in a form that we can input into SCS.  To recover the matrix solution
from the optimal solution returned by SCS, we simply use :math:`S^\star =
\text{mat}(s^\star)`.
