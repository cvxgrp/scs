.. _solver:

Solver
===============
TODO

Primal or dual infeasibility
----------------------------

Cones
-----
The cone K can be any Cartesian product of the following primitive cones.
**Note**:The rows of the data matrix :code:`A` correspond to the cones in :code:`K`. The rows of
:code:`A` must be in the order of the cones given above, i.e., first come the rows that
correspond to the zero/free cones, then those that correspond to the positive
orthants, then SOCs, etc.

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Name
     - Description
     - Entries in :code:`ScsCone`

   * - Zero cone
     - :math:`\{s \mid s = 0 \}`
     - :code:`z`, length of cone
   * - Positive orthant
     - :math:`\{s \mid s \geq 0 \}`
     - :code:`l`, length of cone
   * - Box cone
     - :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^n \mid t l \leq s \leq t u  \}`
     - :code:`bu, bl, bsize`, :math:`u,l` arrays of len :code:`bsize`
   * - Second-order cone
     - :math:`\{(t, s) \in \mathbf{R} \times \mathbf{R}^n\mid \|s\|_2 \leq t  \}`
     - :code:`q, qsize`, array of SOC sizes of len :code:`qsize`
   * - Positive semidefinite cone
     - :math:`\{ S \in \mathbf{R}^{n \times n} \mid \min_i \lambda_i(S) \geq 0, S = S^\top  \}`
     - :code:`s, ssize` see note below
   * - Exponential cone
     - :math:`\{   (x,y,z) \in \mathbf{R}^3 \mid y e^{x/y} \leq z, y>0  \}`
     - :code:`ep`, number of cone triples
   * - Dual exponential cone
     - :math:`\{  (u,v,w)\in \mathbf{R}^3 \mid −u e^{v/u} \leq e w, u<0 \}`
     - :code:`ed`, number of cone triples
   * - Power cone
     - :math:`\{  (x,y,z) \in \mathbf{R}^3 \mid x^p y^{1-p} \geq |z|, x \geq 0, y\geq 0\}`
     - :code:`p, psize` array of powers
   * - Dual power cone
     - :math:`\{ (u,v,w)\in \mathbf{R}^3 \mid (u/p)^p (v/(1-p))^{1-p} \geq |w|, u\geq 0, v\geq 0 \}`
     - :code:`p, psize`

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


Settings
--------

.. list-table::
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Permitted values
     - Default
   * - :code:`normalize`
     - :code:`scs_int`
     - Whether to perform heuristic data rescaling.
     - True/False
     - 1
   * - :code:`init_scale`
     - :code:`scs_float`
     - Initial dual :code:`scale` factor (the important one).
     - :math:`(0, \infty)`
     - 0.1
   * - :code:`adaptive_scaling`
     - :code:`scs_int`
     - Whether to heuristically adapt dual :code:`scale` through the solve.
     - True/False
     - 1
   * - :code:`rho_x`
     - :code:`scs_float`
     - Primal scale factor (less important one).
     - :math:`(0, \infty)`
     - 1e-6
   * - :code:`max_iters`
     - :code:`scs_int`
     - Maximum number of iterations to run.
     - :math:`\mathbf{N}`
     - 1e5
   * - :code:`eps_abs`
     - :code:`scs_float`
     - Absolute feasibility tolerance.
     - :math:`(0, \infty)`
     - 1e-4
   * - :code:`eps_rel`
     - :code:`scs_float`
     - Relative feasibility tolerance.
     - :math:`(0, \infty)`
     - 1e-4
   * - :code:`eps_infeas`
     - :code:`scs_float`
     - Infeasibility tolerance (primal and dual).
     - :math:`(0, \infty)`
     - 1e-7
   * - :code:`alpha`
     - :code:`scs_float`
     - Douglas-Rachford relaxation parameter.
     - :math:`(0, 2)`
     - 1.5
   * - :code:`time_limit_secs`
     - :code:`scs_float`
     - Time limit for solve run in seconds (can be fractional). :code:`0` is interpreted as no limit.
     - :math:`[0, \infty)`
     - 0
   * - :code:`verbose`
     - :code:`scs_int`
     - Whether to print solver output to stdout.
     - True/False
     - 1
   * - :code:`warm_start`
     - :code:`scs_int`
     - Set to True if you initialize the solver with a guess of the solution (see below).
     - True/False
     - 0
   * - :code:`acceleration_lookback`
     - :code:`scs_int`
     - How much memory to use for Anderson acceleration. More memory requires more time to compute but can give more reliable steps. :code:`0` disables it.
     - :math:`\mathbf{N}`
     - 0
   * - :code:`acceleration_interval`
     - :code:`scs_int`
     - Run Anderson acceleration every this number of iterations.
     - :math:`\mathbf{N}`
     - 1
   * - :code:`write_data_filename`
     - :code:`char *`
     - If this is set the problem data is dumped to this filename.
     - Any filename
     - NULL
   * - :code:`log_csv_filename`
     - :code:`char *`
     - If this is set SCS will write csv logs of various quantities through the solver. Doing this makes the solver much slower.
     - Any filename
     - NULL


Termination criteria
--------------------

Acceleration
------------

Warm-starting
-------------

SCS supports warm-starting the solver with a guess of the solution.

Compile flags
-------------
