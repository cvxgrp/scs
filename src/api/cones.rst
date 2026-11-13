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
   * - Complex positive semidefinite cone
     - :math:`\{ s \in \mathbf{R}^{k^2} \mid \text{mat}(s) \succeq 0 \}` (See :ref:`note <sdcone>`)
     - :code:`cs` array of PSD cone lengths with :code:`cssize` elements, each :math:`cs_i = k_i` in description.
     - :math:`\displaystyle \sum_{i=1}^{\text{cssize}} cs_i^2`
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

.. _spectral_cones:

Spectral cones (experimental)
-----------------------------

.. note::

   Spectral cones are **experimental** and must be explicitly enabled at compile
   time by setting :code:`USE_SPECTRAL_CONES=1` (which also requires
   :code:`USE_LAPACK=1`). For example:

   .. code-block:: bash

      make USE_LAPACK=1 USE_SPECTRAL_CONES=1

   or with CMake:

   .. code-block:: bash

      cmake -DUSE_LAPACK=ON -DUSE_SPECTRAL_CONES=ON ..

   The flag also changes the public header: the spectral cone fields of
   :code:`ScsCone` (:code:`d`, :code:`dsize`, :code:`nuc_m`, :code:`nuc_n`,
   :code:`nucsize`, :code:`ell1`, :code:`ell1_size`, :code:`sl_n`,
   :code:`sl_k`, :code:`sl_size`) exist only when :code:`USE_SPECTRAL_CONES`
   is defined, so the struct is larger with it. Any code that includes
   :code:`scs.h` and links against a library built with spectral cones must
   therefore be compiled with :code:`-DUSE_SPECTRAL_CONES` as well. Neither
   build system does this for you: the installed :code:`scs.h` is the same
   file in both cases, and the CMake targets do not pass the definition on to
   the projects that link them. This applies even to programs that use none of
   the spectral cones: the library reads these fields each time a problem is
   set up (in :code:`scs_init`, which :code:`scs` calls), so they must be
   present and set to zero (:code:`SCS_NULL` for the arrays) when unused. A program compiled without the definition passes a struct that is
   too short, and the library reads past the end of it.

   The API and behavior of these cones may change in future releases.
   See :ref:`compile_flags` for more details.

The following cones are available when spectral cone support is enabled.

.. list-table::
   :header-rows: 1

   * - Name
     - Description
     - Entries in :ref:`ScsCone <ScsCone>`
     - Number of rows in :math:`A, b`
   * - Log-determinant cone
     - :math:`\{(t, v, X) \mid -v \log\det(X/v) \leq t,\; v > 0,\; X \succ 0 \}`
     - :code:`d` array of matrix dimensions with :code:`dsize` elements, each :math:`d_i = n_i` where :math:`X \in \mathbf{R}^{n_i \times n_i}` (See :ref:`note <spectral_vec>`).
     - :math:`\displaystyle \sum_{i=1}^{\text{dsize}} \left(\frac{d_i(d_i+1)}{2} + 2\right)`
   * - Nuclear norm cone
     - :math:`\{(t, X) \mid \|X\|_* \leq t \}`
     - :code:`nuc_m`, :code:`nuc_n` arrays of matrix dimensions with :code:`nucsize` elements. :math:`X \in \mathbf{R}^{m_i \times n_i}` with :math:`m_i \geq n_i \geq 1` (See :ref:`note <spectral_vec>`).
     - :math:`\displaystyle \sum_{i=1}^{\text{nucsize}} (m_i n_i + 1)`
   * - :math:`\ell_1` norm cone
     - :math:`\{(t, x) \mid \|x\|_1 \leq t \}`
     - :code:`ell1` array of vector dimensions with :code:`ell1_size` elements.
     - :math:`\displaystyle \sum_{i=1}^{\text{ell1\_size}} (\text{ell1}_i + 1)`
   * - Sum-of-largest-eigenvalues cone
     - :math:`\{(t, X) \mid \sum_{i=1}^{k} \lambda_i(X) \leq t \}`
     - :code:`sl_n`, :code:`sl_k` arrays of matrix dimensions and :math:`k` values with :code:`sl_size` elements, with :math:`0 < k_i < n_i` (See :ref:`note <spectral_vec>`).
     - :math:`\displaystyle \sum_{i=1}^{\text{sl\_size}} \left(\frac{n_i(n_i+1)}{2} + 1\right)`

.. _spectral_vec:

The matrix :math:`X` in these cones must be vectorized as follows. In the
log-determinant and sum-of-largest-eigenvalues cones :math:`X` is symmetric and
is vectorized as for the :ref:`semidefinite cone <sdcone>`: the lower triangle
is stacked column-wise with the off-diagonal entries scaled by :math:`\sqrt{2}`,
which gives a vector of length :math:`n_i(n_i+1)/2`. In the nuclear norm cone
:math:`X` is a general :math:`m_i \times n_i` matrix and its entries are stacked
column by column with no scaling, which gives a vector of length
:math:`m_i n_i`. The scalar entries come first: the rows are stacked as
:code:`[t; v; vec(X)]` for the log-determinant cone and as :code:`[t; vec(X)]`
for the nuclear norm and sum-of-largest-eigenvalues cones.

The nuclear norm cone requires :math:`m_i \geq n_i`. Since
:math:`\|X\|_* = \|X^\top\|_*` a wide matrix can be handled by passing its
transpose. The sum-of-largest-eigenvalues cone requires :math:`0 < k_i < n_i`
(for :math:`k_i = n_i` the constraint is the linear inequality
:math:`\text{Trace}(X) \leq t`). SCS reports a cone validation error at setup
if these conditions do not hold.

The log-determinant cone is the epigraph of the perspective of
:math:`-\log\det`, so with :math:`v = 1` minimizing :math:`t` maximizes
:math:`\log\det X`.


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

The real positive semidefinite cone of matrices is the set

.. math::
   \{S \in \mathbf{R}^{k \times k} \mid  S = S^\top,\  x^\top S x \geq 0 \ \forall x \in \mathbf{R}^k \},

and the complex positive semidefinite cone is the set

.. math::
   \{S \in \mathbf{C}^{k \times k} \mid  S = S^\dagger,\  x^\dagger S x \geq 0 \ \forall x \in \mathbf{C}^k \},

and for short we use :math:`S \succeq 0` to denote membership. SCS
vectorizes these cones in a special way which we detail here.

SCS assumes that the input data corresponding to positive semidefinite cones have been vectorized by scaling the off-diagonal entries by
:math:`\sqrt{2}` and stacking the real parameters of the lower triangle column-wise. For a :math:`k \times k`
real matrix variable (or data matrix) this operation creates a vector of length
:math:`k(k+1)/2`, whereas for a :math:`k \times k` complex matrix it creates a vector of length :math:`k^2`. Scaling by :math:`\sqrt{2}` is required to preserve the inner-product.

**This must be done for the rows of both** :math:`A` **and** :math:`b` **that correspond to semidefinite cones and must be done independently for each semidefinite cone.**

More explicitly, we want to express :math:`\text{Trace}(Y^\dagger S)` as :math:`\text{vec}(Y)^\dagger \text{vec}(S)`,
where the :math:`\text{vec}` operation takes the :math:`k \times k` matrix

.. math::
  S =  \begin{bmatrix}
          S_{11} & S_{12} & \ldots & S_{1k}  \\
          S_{21} & S_{22} & \ldots & S_{2k}  \\
          \vdots & \vdots & \ddots & \vdots  \\
          S_{k1} & S_{k2} & \ldots & S_{kk}  \\
        \end{bmatrix}

and produces a vector consisting of the lower triangular real parameters scaled and re-arranged. In the real case, the result is

.. math::
  \text{vec}(S) = (S_{11}, \sqrt{2} S_{21}, \ldots, \sqrt{2} S_{k1}, S_{22}, \sqrt{2}S_{32}, \dots, S_{k-1,k-1}, \sqrt{2}S_{k,k-1}, S_{kk}) \in \mathbf{R}^{k(k+1)/2},

whereas in the complex case the result is

.. math::
  \text{cvec}(S) ={} & (S_{11}, \sqrt{2} \Re(S_{21}), \sqrt{2} \Im(S_{21}), \ldots, \sqrt{2} \Re(S_{k1}), \sqrt{2} \Im(S_{k1}), S_{22}, \\
  & \sqrt{2}\Re(S_{32}), \sqrt{2}\Im(S_{32}), \dots, S_{k-1,k-1}, \sqrt{2}\Re(S_{k,k-1}), \sqrt{2}\Im(S_{k,k-1}), S_{kk}) \in \mathbf{R}^{k^2}.

To recover the matrix solution this operation must be inverted on the components
of the vectors returned by SCS corresponding to each semidefinite cone. That is, the
off-diagonal entries must be scaled by :math:`1/\sqrt{2}` and the upper triangular
entries are filled in by copying the values of lower triangular entries.
Explicitly, in the real case the inverse operation takes vector :math:`s \in
\mathbf{R}^{k(k+1)/2}` and produces the matrix

.. math::
  \text{mat}(s) =  \begin{bmatrix}
                    s_{1} & s_{2} / \sqrt{2} & \ldots & s_{k} / \sqrt{2}  \\
                    s_{2} / \sqrt{2} & s_{k+1} & \ldots & s_{2k-1} / \sqrt{2}  \\
                    \vdots & \vdots & \ddots & \vdots  \\
                    s_{k} / \sqrt{2} & s_{2k-1} / \sqrt{2} & \ldots & s_{k(k+1) / 2}  \\
                    \end{bmatrix}
  \in \mathbf{R}^{k \times k},

whereas in the complex case the inverse operation takes vector :math:`s \in
\mathbf{R}^{k^2}` and produces the matrix

.. math::
  \text{cmat}(s) =  \begin{bmatrix}
                    s_{1} & (s_{2} - i s_3) / \sqrt{2} & \ldots & (s_{2k-2} - is_{2k-1}) / \sqrt{2}  \\
                    (s_{2} + i s_3) / \sqrt{2} & s_{2k} & \ldots & (s_{4k-5} - is_{4k-4}) / \sqrt{2} \\
                    \vdots & \vdots & \ddots & \vdots  \\
                    (s_{2k-2} + is_{2k-1}) / \sqrt{2} & (s_{4k-5} + is_{4k-4}) / \sqrt{2} & \ldots & s_{k^2}  \\
                    \end{bmatrix}
  \in \mathbf{C}^{k \times k},

So the cone definitions that SCS uses are

.. math::
   \mathcal{S}_\mathbf{R}^k &= \{ \text{vec}(S) \mid S \succeq 0\} = \{s \in \mathbf{R}^{k(k+1)/2} \mid \text{mat}(s) \succeq 0 \}.\\
   \mathcal{S}_\mathbf{C}^k &= \{ \text{cvec}(S) \mid S \succeq 0\} = \{s \in \mathbf{R}^{k^2  } \mid \text{cmat}(s) \succeq 0 \}.

Example
^^^^^^^

For a concrete example in Python see :ref:`py_mat_completion`.
Here we consider the real positive semidefinite cone constraint over
variables :math:`x \in \mathbf{R}^n` and :math:`S \in \mathbf{R}^{k \times k}`

.. math::
    B - \sum_{i=1}^n \mathcal{A}_i x_i = S \succeq 0

where data :math:`B, \mathcal{A}_1, \ldots, \mathcal{A}_n \in \mathbf{R}^{k
\times k}` are symmetric. We can write this in the canonical form over a new
variable :math:`s \in \mathcal{S}_+^k`:

.. math::
  s &= \text{vec}(S)\\
    &= \text{vec}(B - \sum_{i=1}^n \mathcal{A}_i x_i) \\
    &= \text{vec}(B) - \sum_{i=1}^n \text{vec}(\mathcal{A}_i) x_i \\
    &= b - Ax

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
