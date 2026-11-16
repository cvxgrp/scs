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

One solve is treated differently. At the start of every solve, and whenever
the scaling :math:`R` changes, SCS solves the system once for the right-hand
side :math:`(c, -b)` and uses that solution in every iteration until the next
change, so an error in it is repeated in every iteration rather than averaged
out. SCS therefore checks that solve, whatever the linear solver:
it computes the residual from the problem data and, while the residual is
larger than :math:`10^{-12}` relative to the size of the terms
(:math:`10^{-5}` in single precision), it refines the solution by solving for
the residual, for at most three steps. This costs one multiply by the matrix
each time, plus one solve and one multiply per refinement step; with the
indirect solvers a refinement step is another conjugate gradient solve.

The solve in each iteration is checked too, more cheaply. Every 25 iterations
SCS computes the residual of that iteration's solve from the problem data and
compares it with the smaller of SCS's own primal and dual residuals at its last
convergence check (but at least :math:`10^{-12}`, :math:`10^{-5}` in single
precision). If the solve's error is the larger, the iterates cannot become
more accurate than the solve, so SCS refines it with one more solve, for the
residual, and then checks on every iteration until the error is down to a
tenth of the residuals. (Stopping as soon as the error was back below the
residuals made refinement switch off and on again, which slowed some solves
and stalled a few.) Nothing needs to be set for this; what to expect:

- At the default tolerances in double precision it almost never fires (in 9
  of 336 benchmark runs at :math:`10^{-4}` and :math:`10^{-6}`, on 0.01% of
  the solves). The check alone costs one multiply by the matrix every 25
  iterations: about 1% of the linear-system time with the default sparse
  direct solver and about 2% with cuDSS.
- At tight tolerances such as :math:`10^{-8}` or :math:`10^{-9}` the error of a
  direct solve can be as large as the residuals SCS has to reach. That is
  where it fires, and it is what lets such solves converge. While it is on,
  each iteration costs one more multiply and one more solve. On a problem
  that stalls it can stay on for most of the remaining iterations (in the
  benchmarks, SDPs that stalled at :math:`10^{-6}`).
- In single precision the error of the default sparse direct solve can be as
  large as the residuals at the default tolerances, so it can fire there too,
  and on problems that stall it often stays on for most of the solve.
- The indirect solvers already solve to a tolerance tied to these residuals,
  so with them it fires only when a conjugate gradient solve misses that
  tolerance.


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

.. _gpu_reproducibility:

Reproducibility of the GPU solvers
----------------------------------

With either GPU solver, two runs of the same problem on the same machine are
not guaranteed to give bit for bit the same result. The GPU indirect solver
does its sparse matrix-vector products with cuSPARSE's
:code:`CUSPARSE_SPMV_CSR_ALG1`, which cuSPARSE does not guarantee to give the
same bits on every run, and SCS runs cuDSS without its deterministic mode
(:code:`CUDSS_CONFIG_DETERMINISTIC_MODE`), which made LPs and QPs several
times slower to solve on an NVIDIA GB10. SCS has no setting to change either
choice. When two runs differ, their iterates differ slightly. Printed values
that are small near the end of a solve, such as the dual residual and the gap,
can then differ in every digit, and the iteration count can change if a
residual ends close to its tolerance. In a check on an NVIDIA GB10 (CUDA 13.0,
cuDSS 1.0.0), three runs of a random SOCP with each GPU solver all took 150
iterations, and their last iterations differed in the dual residual and the
gap (the gap of the GPU indirect solver by a factor of 13).

.. _new_linear_solver:

Implementing a new linear solver
--------------------------------

In order to implement you own linear system solver, you need to implement the
struct :code:`ScsLinSysWork` that contains the workspace your solver requires,
and implement the functions in :code:`include/linsys.h` as detailed below.
See :code:`linsys` directory for examples.

.. doxygentypedef:: ScsLinSysWork

.. doxygenfile:: linsys.h

