.. _linear_solver:

Linear System Solver
====================

At each iteration :math:`k` SCS solves the following set of linear equations:

.. math::
  \begin{bmatrix}
  R_x + P  &  A^\top \\
  A &  -R_y   \\
  \end{bmatrix}
  \begin{bmatrix}
  x \\
  y
  \end{bmatrix}
  =
  \begin{bmatrix}
  z^k_x \\
  z^k_y
  \end{bmatrix}

for a particular right-hand side :math:`z^k \in \mathbf{R}^{n+m}`. The presence
of the diagonal scaling :math:`R \in \mathbf{R}^{(n+m) \times (n+m)}` matrix is
explained in :ref:`scaling`. The :math:`R_y` term is negated to make the matrix
quasi-definite; we can recover the solution to the original problem from this
modified one. Note that the matrix does not change from iteration to iteration,
which is a major advantage of this approach. Moreover, this system can be solved
approximately so long as the errors satisfy a summability condition, which
permits the use of approximate solvers that can scale to very large problems.
SCS also solves the system once for the right-hand side :math:`(c, -b)` at the
start of every solve and whenever :math:`R` changes; how accurately SCS makes
both kinds of solve is described in :ref:`iterative_refinement`.

.. _iterative_refinement:

Iterative refinement
--------------------

Write :math:`K` for the matrix above. How accurately SCS can solve a problem
is limited by how accurately it solves systems with :math:`K`: its residuals
cannot become much smaller than the error of those solves. Inexact solves are
fine while the residuals are large (the algorithm tolerates errors that
satisfy the summability condition above), but near the requested accuracy an
inexact solve holds the iterates back. A factorization can lose accuracy on a
badly scaled or ill-conditioned :math:`K`, and conjugate gradients stops at a
tolerance. So SCS checks the solves it depends on, by computing their
residuals, and when one is not accurate enough it refines it: it solves with
:math:`K` for the residual and adds the correction to the solution that
residual belongs to. For the :math:`(c, -b)` solve and with the indirect solvers that is an
extra solve; for the direct solvers' per-iteration solves it takes the place
of the plain solve (see below). Two solves are checked, each in its own way:

- the solve for the right-hand side :math:`(c, -b)`, made at the start of a
  solve and whenever :math:`R` changes and used in every iteration until the
  next change, is checked every time it is made against a fixed relative
  threshold (:math:`10^{-12}`, :math:`10^{-5}` in single precision), and
  refined for at most three steps while it is above it (see
  :ref:`g_refinement`);
- each iteration's solve is checked every 25 iterations (every iteration
  while an indirect solver is refining) against SCS's own residuals, and made
  more accurate while its error would hold the iterates back (see
  :ref:`refinement`).

There is nothing to set: both always run, with every linear solver (direct
and indirect, on the CPU and the GPU) and in every interface, and no setting
turns them off. They are deterministic: every decision depends only on the
iterates and the residuals, never on timings or load, so they add no
run-to-run variation of their own. Their costs are given in the two
subsections below.

This is not the per-row refinement of the scaling :math:`R` (see
:ref:`per_row_refinement`), which changes :math:`R` itself, nor the refinement
of the small least-squares solve inside :ref:`acceleration`, nor the
refinement cuDSS does inside its own solves in single precision (see
:ref:`cudss_solver`).

.. _g_refinement:

The solve for :math:`(c, -b)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**When and why.** At the start of every solve, and whenever the scaling
:math:`R`, and so :math:`K`, changes (see :ref:`updating_scale`), SCS solves
:math:`K g = (c, -b)`. Every iteration until the next change uses that
:math:`g`: it enters the quadratic that gives :math:`\tau` (the root-plus
function in :ref:`scaling`, where it is written :math:`r`), and :math:`\tau g`
is subtracted from that iteration's solve. An error in an iteration's own
solve affects one iteration; an error in :math:`g` is repeated in every
iteration rather than averaged out, and can limit the accuracy SCS reaches. So
:math:`g` is checked every time it is computed, against a fixed small
threshold, rather than against SCS's current residuals as each iteration's
solve is.

**The check.** After the solve, SCS computes the residual :math:`r = (c, -b) -
K g` with its own products from the problem data, also with a linear solver
that computes the per-iteration products itself (see "Who computes the
multiplies" in :ref:`refinement`), and divides :math:`\|r\|_\infty` by the
size of the terms whose rounding errors it is compared with: the largest
:math:`\infty`-norm among :math:`(c, -b)`, :math:`R g`, :math:`A^\top g_y`,
:math:`P g_x + A^\top g_y` and :math:`A g_x`. If this relative residual is
above :math:`10^{-12}` (:math:`10^{-5}` in single precision, which cannot
reach :math:`10^{-12}`; ``G_SOLVE_TOL`` in ``include/glbopts.h``), SCS refines
:math:`g`.

**The refinement.** Each step solves :math:`K d = r` and computes the residual
of :math:`g + d`. SCS keeps :math:`g + d` unless its residual is larger (or
not finite), in which case it keeps :math:`g` and stops. It also stops once
the residual is at most the threshold, after a step that does not at least
halve the residual, or after three steps (``G_REFINE_MAX_STEPS``). So the
refined :math:`g` never has a larger residual than the first solve.

**Cost.** One multiply by :math:`K` each time :math:`g` is solved for, paid
even when the first solve is already accurate enough, plus one solve and one
multiply per refinement step. With a direct solver each solve reuses the
existing factorization, so this is cheap next to factorizing. With the
indirect solvers the solve for :math:`g` is a conjugate gradient solve
started from zero, to the solvers' best accuracy (:code:`CG_BEST_TOL`).
Conjugate gradients stops on the absolute residual of its reduced system,
which on badly scaled data can leave a much larger relative residual here, so
:math:`g` is checked and refined with them too. There each refinement step is
another conjugate gradient solve from zero, which can take up to about as
many conjugate gradient iterations as the solve for :math:`g` itself.

.. _refinement:

Accuracy of each iteration's solve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The right-hand side of each iteration's solve changes every iteration, and the
solve's error only matters once it is comparable to SCS's residuals. So SCS
checks it periodically and makes the solves more accurate only while their
error would hold the iterates back.

**The check.** Every 25 iterations SCS computes the residual of that
iteration's solve, :math:`\|z^k - K u\|_\infty` with :math:`K` the matrix
above, and compares it with *tol*, the smaller of SCS's own primal and dual
residuals at its last convergence check (but at least :math:`10^{-12}`,
:math:`10^{-5}` in single precision). While refinement is off, the check costs
one multiply by :math:`K` every 25 iterations.

**Direct solvers** (sparse direct, dense, Apple Accelerate, MKL Pardiso,
cuDSS). When the residual of a check is above :math:`0.01\,` *tol*, each
following iteration solves in *correction form*: from the previous
iteration's solution :math:`u_0` it solves :math:`K d = z^k - K u_0` and takes
:math:`u = u_0 + d`. The error of a solve with a factorization is about
proportional to the size of the right-hand side it is given. The iterates
move little from one iteration to the next, so :math:`z^k - K u_0` is
:math:`10^{-3}` to :math:`10^{-7}` of :math:`z^k`, and :math:`u` is that much
more accurate than a plain solve, for one more multiply by :math:`K` per
iteration and no extra solve. The errors do not build up from one iteration to
the next: each correction-form solve shrinks the error it inherits by the same
factor as it shrinks the right-hand side. When the previous solution is no
better a guess than zero (:math:`\|z^k - K u_0\| \ge \|z^k\|`, which
happens after :math:`R` changes) that iteration solves plainly. At each check
while it is on, SCS also makes a plain solve of the same right-hand side, only
to measure its error, and turns the correction form off once that error is at
most :math:`0.001\,` *tol*; the iterate stays the correction-form solution,
unless the correction-form solution's own residual (one more multiply by
:math:`K` per check) is more than twice the plain solve's and above rounding
level (:math:`10^{-14}\,\|z^k\|_\infty` in double precision,
:math:`10^{-5}\,\|z^k\|_\infty` in single). That safeguard is for a
factorization so inaccurate that the correction form does not reduce the
error; normally the correction-form residual is near rounding level and the
safeguard changes nothing.

Why these choices: what slows SCS down is switching between accurate and
inexact solves, more than the inexact solves themselves. On QFORPLAN at
:math:`10^{-6}` with Apple Accelerate, the previous rule (refine only once the
error exceeds *tol*, with one more solve per iteration, until it falls below
:math:`0.1\,` *tol*) took 1,300 to 140,900 iterations depending on the
starting point, and the correction form with the plain solve put into the
iterate at every check (one inexact solve every 25 accurate ones) took up to
31 times the iterations of the rule here, which took 700 to 975. Turning on at
:math:`0.01\,` *tol* and off at :math:`0.001\,` *tol* starts the accurate
solves before they are needed and keeps them while the residuals fall.
Thresholds ten times lower were on for 12.6% of the iterations at the default
tolerances instead of 1.3%, and for most of the solve on more problems, with
similar results at tight tolerances; ten times higher ones took more
iterations on the badly scaled problems.

**Indirect solvers** (CPU and GPU conjugate gradients). Conjugate gradients
already solves to a tolerance tied to SCS's residuals, and its error is set by
that tolerance rather than by the size of the right-hand side, so the
correction form would not make it more accurate. When a check finds the
residual above *tol*, that iteration's solve and each following one is
followed by one refinement step (one more conjugate gradient solve, for the
residual). While this is on, every iteration is a check (one multiply by
:math:`K`, of the solve before its refinement step), and it turns off at the
first that finds the residual at most :math:`0.1\,` *tol*. This fires only
when a conjugate gradient solve misses its tolerance.

**What to expect** (direct solvers):

- At the default tolerances it rarely turns on: on 112 LPs, QPs and SDPs at
  :math:`10^{-4}` and :math:`10^{-6}` from 3 starting points each (QDLDL), for
  1.3% of the iterations, with the same iterations in 318 of the 336 runs (9
  fewer, 9 more) and about 1% more time, within the run-to-run noise.
- At tight tolerances it is what lets direct solves converge, and it takes
  fewer iterations than refining once the error is already too large: on 500
  random LPs and SOCPs at :math:`10^{-9}` with QDLDL, 2.6% fewer iterations and
  2.3% less time; with cuDSS on 200 of them, 7% fewer iterations and 9% less
  solve time (14% less without Anderson acceleration); with Apple Accelerate,
  whose solves can be inaccurate on badly scaled problems, on 30 LPs, QPs and
  SDPs at :math:`10^{-4}` and :math:`10^{-6}` from 3 starting points each,
  0.95 times the time in geometric mean of the per-run ratios, and 0.89 times
  the shifted geometric mean of the times (shift 10 s), which weighs the long
  runs more (QFORPLAN at :math:`10^{-6}`: 700 to 975 iterations, where the
  previous rule took 1,300 to 140,900 depending on the starting point). On 45 LPs and QPs at :math:`10^{-8}` and
  :math:`10^{-9}` with QDLDL the results were the same within noise. MKL
  Pardiso's solves are accurate enough that it turns on less often: on 33 LPs
  and QPs at :math:`10^{-9}`, in 15 of them, for 3.5% fewer iterations and the
  same time.
- While it is on, an iteration's linear-system time is about 1.1 to 1.5 times
  that of a plain one: the extra multiply by :math:`K` costs most, relative to
  the solve, on small and very sparse systems (1.4 with QDLDL on small LPs,
  1.1 with Accelerate on larger QPs).
- On a problem that stalls it can stay on for most of the remaining
  iterations: on 26 SDPs and QPs that stall at :math:`10^{-4}` to
  :math:`10^{-6}` (78 runs), it was on for 19% of the iterations and for most
  of the solve in 11 runs, which solved one run more and took 6% more time.
  In single precision it tends to stay on once it turns on, because the
  threshold to turn it off is near or below what a single-precision solve can
  reach.

**Who computes the multiplies.** By default SCS computes the multiplies by
:math:`K` itself, from the problem data, on one thread. A linear solver can
compute them with its own copy of :math:`K` instead (see
:ref:`new_linear_solver`): MKL Pardiso uses MKL's threaded sparse matrix-vector
product, and cuDSS computes them on the GPU, where the vectors stay. The
decisions above are the same either way; only the cost of the multiply and its
rounding change.

What the solver's own multiplies buy is a cheaper iteration while the
correction form is on. With cuDSS on 10 large LPs, QPs and SOCPs at
:math:`10^{-8}` (NVIDIA GB10), such an iteration cost 1.19 times a plain one
with the multiplies on the GPU, against 1.48 with SCS's on the CPU (on BOYD1
1.61 ms against 2.63, a plain iteration taking 1.33). With MKL Pardiso on 33
LPs and QPs at :math:`10^{-9}` (4 threads, AMD Zen 3) it cost 1.32 times a plain
one against 1.45. Whole solves changed by a few percent at most, within the
run-to-run variation of the iteration counts, because the correction form is
on for a small part of most solves.


Available linear solvers
------------------------

Each of the below linear solvers is included in their own binary. If linking
against SCS directly, then to switch between them you must compile and link
against the right binary. If calling SCS via one of the interfaces then you can
choose between the different linear solvers using the appropriate settings.

.. _direct:

Sparse direct
^^^^^^^^^^^^^

The direct method is the default linear system solver in SCS and factorizes the
above matrix using a sparse permuted LDL factorization. Then it solves the
linear system at each iteration using the cached factors.  Since the linear
system is quasidefinite we have strong existence guarantees about the
factorization.  It relies on the external (but included) `AMD
<https://github.com/DrTimothyAldenDavis/SuiteSparse>`_ and `QDLDL
<https://github.com/oxfordcontrol/qdldl>`_ packages.

.. _dense:

Dense direct
^^^^^^^^^^^^

The dense direct method reduces the KKT system to the smaller Gram matrix
:math:`G = R_x + P + A^\top R_y^{-1} A` (size :math:`n \times n`) and
factorizes it using the LAPACK ``dpotrf`` (Cholesky) routine. The
:math:`A^\top R_y^{-1} A` product is computed via ``dsyrk`` and subsequent
solves use ``dpotrs``. When the diagonal :math:`R` changes, the Gram matrix is
re-formed and re-factorized.

This backend is best suited for small to medium-sized problems where the
constraint matrix :math:`A` is dense. For such problems, dense BLAS/LAPACK
routines can outperform sparse solvers due to highly optimized memory access
patterns and lower overhead. It requires LAPACK (``USE_LAPACK=1``).

.. _mkl:

MKL Pardiso
^^^^^^^^^^^
The `Intel MKL Pardiso solver
<https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface.html>`_
is a shared-memory multiprocessing parallel direct sparse solver which is
included Intel oneAPI MKL library. This offers an alternative to the single
threaded AMD / QDLDL libraries which come bundled with SCS. Pardiso tends to be
faster than AMD / QDLDL, especially for larger problems. If MKL is installed on
your system then it is generally worth using MKL for both the blas / lapack
usage as well as the linear system solve.
Intel MKL is now available for
`free and without restrictions for everyone <https://www.intel.com/content/www/us/en/developer/articles/news/free-ipsxe-tools-and-libraries.html>`_,
though it only offers limited support for non-Intel CPUs. The x86-64 manylinux
Python wheels include it, statically linked (see :ref:`python_install`).

**Threads.** With a threaded MKL, Pardiso runs on MKL's threads
(``MKL_NUM_THREADS``, else ``OMP_NUM_THREADS``): the reordering, the
factorization and the solve done in every iteration, which SCS asks Pardiso to
run in parallel. Starting and synchronizing
the threads costs a few tens of microseconds per call, which on a small linear system
is more than the work. So SCS runs Pardiso on one thread when the factor
Pardiso computes has fewer than 50,000 nonzeros (``PARDISO_MIN_THREADED_NNZ``
in ``linsys/mkl/direct/private.c``). The choice is made once, at setup, from
the problem data and the ordering Pardiso computes; all the Pardiso calls of a
solve then use the same thread count, and MKL's BLAS and LAPACK keep their
threads. On 4 threads, small problems iterated up to about 2x faster this
way with ``MKL_DYNAMIC=FALSE``, and up to about 1.3x with MKL's default
(dynamic mode, which already runs the smallest systems on one thread by
itself); larger ones were unaffected. The break-even point depends on the
problem and the machine; in our measurements it was between about 30,000 and
60,000 factor nonzeros. A sequential MKL, such as the one in the Python
wheels, always runs on one thread, so there this changes nothing. On small
problems a sequential MKL is still faster than a threaded one, because a
threaded MKL costs more per call even on one thread. On several threads
Pardiso may round differently from run to run, so two runs of the same problem
can differ slightly in their results and iteration counts.

**Multiplies for the refinement.** MKL Pardiso computes the multiplies by
:math:`K` of the per-iteration refinement (see :ref:`refinement`) with MKL's
sparse BLAS, on the same upper triangle Pardiso factorizes, analysed once at
setup (``mkl_sparse_optimize``). They run on MKL's threads, or on one thread
when Pardiso does (the rule above), so that a solve never mixes thread counts.
Measured per multiply against SCS's own single-threaded product (AMD Zen 3, 4
threads, 46 LPs and QPs): 0.2 to 0.5 times the time where it runs on MKL's
threads (0.76 against 2.39 ms on co-100), and 0.45 to 0.93 times where it runs
on one thread; forcing 4 threads on the smallest systems made it about twice
as slow as one thread. In the solves, an iteration with the correction form
on cost 1.32 times a plain one, against 1.45 with SCS's product (33 LPs and QPs
at :math:`10^{-9}`). The products add no run-to-run variation of their own:
with the same number of threads they gave the same bits in repeated runs (a
different thread count can change their last bits). The products also cost
setup time that has not been measured: MKL analyses :math:`K`
(``mkl_sparse_optimize``) at setup and again at every update of :math:`R`,
whether or not the correction form ever turns on, since the check every 25
iterations uses the product too. In a build with ``DLONG`` but 32-bit MKL
integers, a :math:`K` too large for those integers (more than :math:`2^{31}`
stored nonzeros) gets SCS's own products instead, as before this change;
``BLAS64=1`` (Meson's ``use_blas64``) provides MKL's at any size. If they
cannot be set up otherwise (out of memory), setup fails with a message.

.. _apple_accelerate:

Apple Accelerate
^^^^^^^^^^^^^^^^

The Apple `Accelerate
<https://developer.apple.com/documentation/accelerate>`_ framework provides
a sparse LDL\ :sup:`T` factorization that is optimized for Apple hardware
(including Apple Silicon). This backend uses the unpivoted LDL\ :sup:`T`
solver from the Accelerate Sparse Solvers API, which is well-suited for the
quasi-definite KKT systems that SCS produces. Accelerate performs its own
fill-reducing ordering internally so no external AMD package is needed.

This backend is macOS-only and is available without installing any additional
libraries since the Accelerate framework ships with Xcode / the macOS SDK.
It does not support 64-bit integer indexing (``DLONG``).

To build with Make::

    make accelerate

To build with CMake::

    cmake -DUSE_APPLE_ACCELERATE=ON ..

In Python, the Accelerate backend is included automatically on macOS
(``pip install scs`` is sufficient). On macOS,
``linear_solver=scs.LinearSolver.AUTO`` selects the bundled QDLDL; select
Accelerate explicitly with::

    solver = scs.SCS(data, cone, linear_solver=scs.LinearSolver.ACCELERATE)

.. _indirect:

Sparse indirect
^^^^^^^^^^^^^^^

The indirect method solves the above linear system approximately with a
'matrix-free' method. To do this it first reduces the system to solving

.. math::

  (R_x + P + A^\top R_y^{-1} A) x & = z^k_x + A^\top R_y^{-1} z^k_y \\
                            y & = R_y^{-1}(A x - z^k_y).

then solves the positive definite system using using `conjugate gradients
<https://en.wikipedia.org/wiki/Conjugate_gradient_method>`_.  Each iteration of
CG requires one multiply each of sparse matrices :math:`P, A, A^\top`.  The
system is solved up to some tolerance, which is tuned to ensure that the
overall algorithm converges: each solve targets a small fraction
(:code:`CG_TOL_FACTOR`, default :math:`0.03`) of the smaller of the current
primal and dual residuals, additionally capped by a term decaying with
iteration :math:`k` like :math:`O(1/k^\gamma)` where :math:`\gamma` is
determined by the constant :code:`CG_RATE` (defaults to :math:`2`), and
floored at :code:`CG_BEST_TOL`. Each solve is warm-started from the previous
solution. The solver also harvests approximate eigenvectors of the reduced
matrix from the CG iterations themselves (eigCG; Stathopoulos and Orginos
2010) and uses up to :code:`DEFLATE_VECTORS` (default 30) of them to deflate
the smallest eigenvalues from subsequent solves, which substantially reduces
the CG iterations needed on ill-conditioned problems. Applying the deflation
in a solve costs no matrix multiplies. Each harvest costs the multiplies of
one CG iteration per vector.

Deflation needs LAPACK and double precision. It is switched off when SCS is
built without LAPACK (:code:`USE_LAPACK=0`) or in single precision
(:code:`SFLOAT=1`), and only the CPU indirect solver has it. Single precision
builds also use :math:`0.2` for :code:`CG_TOL_FACTOR` and :math:`1.5` for
:code:`CG_RATE`.

Deflation keeps :math:`2n` numbers per vector in memory (:math:`60n` by
default), and a harvest takes a further :math:`n` numbers per window entry and
:math:`2n` per vector while it runs (:math:`150n` by default). Two environment
variables change this for a run without recompiling. Both are read once, when
the solver is set up (:code:`scs_init`):

* :code:`SCS_DEFLATE` sets the number of deflation vectors in place of
  :code:`DEFLATE_VECTORS`. :code:`SCS_DEFLATE=0` turns deflation off.
* :code:`SCS_EIGCG_WIN` sets the size of the eigCG window, the number of CG
  vectors that a harvest keeps in memory (default 60, :code:`EIGCG_WINDOW` in
  :code:`linsys/cpu/indirect/private.c`). The window is never smaller than
  three times the number of deflation vectors, so with the default 30 vectors
  it is 90 and a smaller :code:`SCS_EIGCG_WIN` has no effect unless
  :code:`SCS_DEFLATE` is lowered as well.

The indirect method has the advantage of not requiring an expensive
factorization but typically is slower on a per-iteration basis. In most cases
the factorization is relatively cheap so the direct method is the default,
however for very large problems the indirect solver can be faster.

.. _gpu_indirect:

Sparse GPU indirect method
^^^^^^^^^^^^^^^^^^^^^^^^^^

The above linear solvers all run on CPU. We also have support for a GPU version
of the indirect solver, where the matrix multiplies are all performed on the
GPU. It runs conjugate gradients with the same diagonal preconditioner,
tolerance schedule and warm start as the CPU indirect solver, but it does not
harvest eigenvectors or deflate, so :code:`DEFLATE_VECTORS` and the
:code:`SCS_DEFLATE` and :code:`SCS_EIGCG_WIN` environment variables have no
effect on it. Usually the direct GPU solver will be faster than this solver.

.. _cudss_solver:

Sparse GPU direct method
^^^^^^^^^^^^^^^^^^^^^^^^^

This is a linear solver that uses the `cuDSS
<https://developer.nvidia.com/cudss>`_ library to solve the linear system on
the GPU. It is similar to the direct solver, but uses the calls to cuDSS
library to perform the analysis, numerical (re-)factorization and subsequent
solves. According to its documentation reordering (a major part of the analysis
phase) is executed on the host while symbolic factorization (another part of
the analysis phase), numerical factorization, and solve are executed on the
GPU.

SCS runs cuDSS without pivoting, matching or iterative refinement, so it
computes a plain LDL' factorization, as the CPU direct solver does. The linear
system is quasi-definite, so this factorization exists under any ordering, and
each iteration costs one solve. From version 1.0, cuDSS by default applies a
matching to a symmetric indefinite matrix. With its default pivoting that
makes the factorization up to four times larger (2.4 times larger on average
on the LPs and QPs measured), and without pivoting it made the solves
inaccurate enough for SCS to stall on some LPs. Older versions of cuDSS do not
apply a matching unless asked to, so SCS factorizes the same way with every
cuDSS version it supports. The tradeoff is that the accuracy of a solve rests
on the conditioning of the system: the worst relative residual of a solve in
the test suite was about 6e-8, against 5e-11 with cuDSS's default pivoting and
matching. Single-precision builds keep one step of iterative refinement per
solve, which about doubles the cost of a solve, because without it one problem
in the test suite does not converge.

**Multiplies for the refinement.** When :math:`K` has at least 60,000 nonzeros
(``CUDSS_MIN_PRODUCT_NNZ`` in ``linsys/cudss/direct/private.c``), the cuDSS
solver computes the multiplies by :math:`K` of the per-iteration refinement
(see :ref:`refinement`) on the GPU: one cuSPARSE sparse matrix-vector product
and one cuBLAS reduction, with the guess and the residual kept on the GPU, so
an iteration in correction form copies no more vectors between the host and
the GPU than a plain solve. cuSPARSE has no symmetric product, so :math:`K`
is stored a second time with both triangles, about twice the nonzeros of the
copy cuDSS factorizes. A call costs about 0.04 ms whatever the size, so smaller
systems use SCS's product on the CPU, which was faster there. On an NVIDIA
GB10, per multiply, the GPU took 0.12 to 0.37 times the CPU's time above
:math:`3.7 \cdot 10^5` nonzeros (0.18 against 1.27 ms on BOYD1).
In the solves, on 10 large LPs, QPs and SOCPs at :math:`10^{-8}`, an
iteration with the correction form on cost 1.19 times a plain one, against
1.48 with the multiplies on the CPU, and the time per iteration fell by 4%.
Setting them up (building :math:`K` with both triangles, copying it to the
GPU, creating the cuSPARSE and cuBLAS handles) is done at setup whether or not
the correction form ever turns on, since the check every 25 iterations uses
the multiply too. On the GB10 it took 3, 35 and 340 ms for random :math:`K`
with :math:`10^5`, :math:`10^6` and :math:`8 \cdot 10^6` nonzeros, 5 to 13%
of the time cuDSS's own setup (analysis and factorization) took; the two
handles took 0.4 to 0.8 ms of that. If they cannot be set up (out of GPU
memory, or a :math:`K` too large for 32-bit indices in a build without
``DLONG``), setup fails with a message.

.. _gpu_reproducibility:

Reproducibility of the GPU solvers
----------------------------------

With either GPU solver, two runs of the same problem on the same machine are
not guaranteed to give bit for bit the same result. The GPU indirect solver
does its sparse matrix-vector products with cuSPARSE's
:code:`CUSPARSE_SPMV_CSR_ALG1`, which cuSPARSE does not guarantee to give the
same bits on every run, and SCS runs cuDSS without its deterministic mode
(:code:`CUDSS_CONFIG_DETERMINISTIC_MODE`), which made LPs and QPs several times
slower to solve on an NVIDIA GB10. SCS has no setting to change either choice.
The multiplies of the refinement on the GPU (cuDSS, see above) use
:code:`CUSPARSE_SPMV_CSR_ALG2`, which gives the same result every run, so they
add no variation of their own. When two runs differ, their iterates differ
slightly. Printed values that are small near the end of a solve, such as the
dual residual and the gap, can then differ in every digit, and the iteration
count can change if a residual ends close to its tolerance. In a check on an
NVIDIA GB10 (CUDA 13.0, cuDSS 1.0.0), three runs of a random SOCP with each GPU
solver all took 150 iterations, and their last iterations differed in the dual
residual and the gap (the gap of the GPU indirect solver by a factor of 13).

.. _new_linear_solver:

Implementing a new linear solver
--------------------------------

In order to implement you own linear system solver, you need to implement the
struct :code:`ScsLinSysWork` that contains the workspace your solver requires,
and implement the functions in :code:`include/linsys.h` as detailed below.
See :code:`linsys` directory for examples.

**Products with the matrix (optional).** The per-iteration refinement of the
direct solvers (see :ref:`refinement`) multiplies by the matrix :math:`K`:
every 25 iterations to check a solve, and once per iteration while the
correction form is on. By default SCS computes these products itself, from the
problem data, on one thread. A solver that holds :math:`K` anyway can compute
them faster with its own copy, by implementing three functions:

- :code:`scs_lin_sys_residual(w, b, x)` computes :math:`r = b - K x` with the
  solver's current :math:`K` (including the latest :math:`R` from
  :code:`scs_update_lin_sys_diag_r`), keeps :math:`r`, and returns
  :math:`\|r\|_\infty`. If :code:`x` is :code:`SCS_NULL`, :math:`x` is the
  vector stored by the last :code:`scs_lin_sys_keep`.
- :code:`scs_lin_sys_keep(w, x)` stores :math:`x` as the guess :math:`x_0`.
  SCS calls it with the solution of the latest solve (a
  :code:`scs_solve_lin_sys` or :code:`scs_lin_sys_correct` call), unchanged,
  so a GPU solver can copy the solution it still has on the device and not
  read :code:`x` at all. Between that solve and this call SCS may call
  :code:`scs_lin_sys_residual` (when it turns the correction form on, it
  checks the solve with :code:`scs_lin_sys_residual(w, b, x)` first), but no
  other solve, so the residual must not overwrite the solver's copy of the
  solution.
- :code:`scs_lin_sys_correct(w, x)` returns :math:`x = x_0 + K^{-1} r` for the
  :math:`r` of the latest :code:`scs_lin_sys_residual(w, b, SCS_NULL)`.

Every solver must define the three functions so that it links (an
out-of-tree solver written for an earlier SCS needs the three stubs below).
The contract:

- :code:`scs_lin_sys_residual` returns :math:`\|r\|_\infty \ge 0`,
  :code:`SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED` if the workspace does not compute
  the products, or :code:`SCS_LIN_SYS_PRODUCTS_FAILED` on an error.
- Whether a workspace provides the products is fixed when
  :code:`scs_init_lin_sys_work` returns, by a rule stated up front (cuDSS: the
  number of nonzeros of :math:`K`; MKL Pardiso: always, except a :math:`K` too
  large for 32-bit MKL integers in a ``DLONG`` build without ``BLAS64``), never
  by whether memory or a library happened to be available. A solver that should provide
  them but cannot set them up makes :code:`scs_init_lin_sys_work` fail with a
  message saying why, rather than returning "not provided".
- SCS calls :code:`scs_lin_sys_residual` once at setup, with zero vectors:
  "not provided" makes SCS compute every product itself (and never call the
  other two functions); a failure makes SCS's setup fail; a norm means the
  solver provides them. During a solve any negative return is a failure of
  the linear solver.

A solver without the products defines::

    scs_float scs_lin_sys_residual(ScsLinSysWork *w, const scs_float *b,
                                   const scs_float *x) {
      return SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED;
    }
    void scs_lin_sys_keep(ScsLinSysWork *w, const scs_float *x) {}
    scs_int scs_lin_sys_correct(ScsLinSysWork *w, scs_float *x) { return -1; }

SCS makes every decision (when to check, the thresholds, when to use the
guess); the solver only computes, so the refinement behaves the same with
every solver, up to the rounding of the product.
:code:`test_lin_sys_products` in the test suite checks a solver's products
against SCS's own. The MKL Pardiso solver (:code:`linsys/mkl/direct/private.c`,
with MKL's sparse BLAS) and the cuDSS solver
(:code:`linsys/cudss/direct/private.c`, on the GPU with cuSPARSE and cuBLAS)
are examples; the other solvers return "not provided".

.. doxygentypedef:: ScsLinSysWork

.. doxygenfile:: linsys.h

