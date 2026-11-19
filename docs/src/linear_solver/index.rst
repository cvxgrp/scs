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
residual belongs to. For the :math:`(c, -b)` solve and with the indirect
solvers that is an extra solve; for the direct solvers' per-iteration solves
it takes the place of the plain solve (see below). Two solves are checked,
each in its own way:

- the solve for the right-hand side :math:`(c, -b)`, made at the start of a
  solve and whenever :math:`R` changes and used in every iteration until the
  next change, is solved to about machine precision every time it is made,
  with at most five refinement steps (see :ref:`g_refinement`);
- each iteration's solve is checked every 25 iterations against SCS's own
  residuals, and made more accurate while its error would hold the iterates
  back (see :ref:`refinement`).

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
:ref:`cudss_solver`), nor the two refinement steps MKL Pardiso takes inside a
solve when its factorization had to perturb a pivot (``iparm[7] = 0`` with
``iparm[9] = 13``, see :ref:`mkl`).

SCS's linear solve takes a warm start, a tolerance and a number of refinement
steps (see :ref:`new_linear_solver`). The per-iteration refinement uses it
for both kinds of solver and decides only when: with the direct solvers it
passes the previous iteration's solution as the warm start, and with the
indirect solvers it asks for refinement steps. The solve for :math:`(c, -b)`
passes its first solution as the warm start, a tolerance at about machine
precision and at most five steps in all, with every solver (see below).

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
:math:`g` is solved to about machine precision every time it is computed,
rather than as accurately as SCS's current residuals need, as each
iteration's solve is.

**The tolerance.** A plain solve (with conjugate gradients: from zero, to
the solvers' best accuracy, :code:`CG_BEST_TOL`) gives :math:`g_0` and its
residual :math:`r = (c, -b) - K g_0`. Computing :math:`r` rounds each row by
a few machine epsilons :math:`\epsilon` of the size of the terms it is
computed from, so SCS asks for

.. math::

  \|(c, -b) - K g\|_\infty \le 4 \epsilon \, \big\| |(c, -b)| + |K|\,|g_0| \big\|_\infty

(absolute values entry by entry; :code:`G_TOL_EPS` in
``include/glbopts.h``), with :math:`\epsilon` the machine epsilon of the
precision SCS is built in. The size is computed row by row, from the
entries each row actually multiplies, rather than as
:math:`\|(c, -b)\|_\infty + \|K\|_\infty \|g_0\|_\infty` (the normwise
backward error), which pairs the largest row of :math:`K` with the largest
entry of :math:`g_0` even where they do not meet: on badly scaled problems
that was up to :math:`10^9` times larger, and would have accepted residuals
far above rounding. It is one number, the largest row's, because the
tolerance of the solve bounds the largest entry of the residual; so a row
whose terms are much smaller than the largest row's is held to the largest
row's rounding, not to its own.

**The refinement.** If :math:`g_0` misses the tolerance, the solve is
repeated with :math:`g_0` as its warm start, to the tolerance, with at most
five steps in all (:code:`G_REFINE_MAX_STEPS`). For a direct solver the
first step is the correction form from :math:`g_0`; for conjugate gradients
it is a solve started from :math:`g_0`. Each further step computes the
residual with :math:`K` and adds the solution of :math:`K d = r` (by
conjugate gradients from zero to the tolerance, with them). A step is kept
only if it lowers the residual, the steps stop after one that does not at
least halve it, and SCS keeps :math:`g_0` unless the result is better, so
:math:`g` never has a larger residual than the first solve. Conjugate
gradients solves a reduced system and computes :math:`g_y` from
:math:`g_x`, so on badly scaled data its own rounding can be far above that
of :math:`K`'s residual: the solve started from :math:`g_0` is limited by
that rounding, but each further step solves only for the correction, so it
is not.

**If it fails.** The solve should never fail. If the linear solver returns
an error, :math:`g` is not finite, or the residual of the refined :math:`g`
is more than half the size of its terms (:code:`G_FAIL_REL_RES`: not even
the leading digit right; with working solvers it measured at most about
:math:`10^{-2}` of the size, and that only where the solver could do no
better, such as conjugate gradients at its iteration limit), SCS stops with
status :code:`SCS_FAILED` (see :ref:`exit_flags`) and prints what failed,
that this indicates a bug in SCS or its linear solver, and a request to open
an issue at https://github.com/cvxgrp/scs/issues with the problem data (the
``write_data_filename`` setting saves it). If :math:`b` or :math:`c`
contains a NaN or infinite entry, SCS stops with :code:`SCS_FAILED` and says
so instead.

**Cost.** One multiply by :math:`K` each time :math:`g` is solved for (the
residual of :math:`g_0`) and one pass over :math:`A` and :math:`P` for the
size of the terms, paid even when the first solve is already accurate
enough; when it is not, one solve and one multiply per step, and with a
direct solver one multiply more (the correction form computes the residual
of :math:`g_0` again). With a direct solver each solve reuses the existing
factorization, so this is cheap next to factorizing. With the indirect
solvers each step costs the conjugate gradient iterations that reach the
tolerance. On random
SOCPs the refinement took up to about 30% of the first solve's conjugate
gradient iterations in double precision, and up to about two thirds in
single precision. Where conjugate gradients cannot reach the tolerance
within its iteration limit (:math:`10 n`), each step runs to that limit, so
on such problems the refinement can cost up to about five times the first
solve.

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
residuals at its last convergence check. While refinement is off, the check
costs one multiply by :math:`K` every 25 iterations.

**Rounding level.** In double precision refinement is never on while that
residual is at rounding level: at most four machine epsilons times
:math:`\|\,|z^k| + |K|\,|u|\,\|_\infty`, the size of the terms the residual
is computed from, row by row. A residual that small is mostly the rounding of
its own computation, and refinement cannot lower it. This matters where *tol*
goes to zero, on an infeasibility or unboundedness certificate: before 3.4.0
the floor was absolute (*tol* at least :math:`10^{-12}`), below rounding
level wherever those terms are large, so there refinement turned on and
stayed on to the end. The size is computed only at a check where it can
change the outcome (one pass over :math:`A` and :math:`P`, the cost of a
multiply). Single precision has no such test: there refinement also lowers an
error of the solution that the residual does not show (with conjugate
gradients and the dense solver, which solve a smaller system), and turning
it off at rounding level lost problems at :math:`10^{-6}`. With conjugate
gradients refinement turns off, in both precisions, once the error is at
most the solvers' best accuracy (:math:`10^{-12}`, :math:`10^{-6}` in single
precision), below which they are never asked to go.

**Direct solvers** (sparse direct, dense, Apple Accelerate, MKL Pardiso,
cuDSS). When the residual of a check is above :math:`0.01\,` *tol*, each
following iteration solves in *correction form*: SCS keeps the previous
iteration's solution :math:`u_0` and passes it to the linear solve as its
warm start, which solves :math:`K d = z^k - K u_0` and takes
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
correction form would not make it more accurate. What goes wrong is that a
solve misses its tolerance: it stops at its iteration limit, the residual it
tracks drifts from the true one (in single precision), or :math:`y`, which it
computes from :math:`x`, carries rounding that its smaller system does not
see. When a check finds the residual above *tol*, each following iteration's
solve also refines itself: it computes the residual of the whole system (one
multiply by :math:`K`) and, while that is above the conjugate gradient
tolerance, solves for it with conjugate gradients from zero, to that
tolerance, and adds the correction, keeping each such step only if it lowers
the residual, and at most two. So refinement costs one multiply per
iteration, and conjugate gradient iterations only where a solve missed its
tolerance. At each check while it is on, SCS also makes a plain conjugate
gradient solve of the same right-hand side, only to measure its error, and
turns refinement off once that error is at most :math:`0.1\,` *tol*, as for
the direct solvers. The steps solve for the residual of the whole system:
restarting conjugate gradients on the same right-hand side from the
solution, the same in exact arithmetic, works with the residual of the
smaller system and cannot see all of the rounding error: in single precision
it solved 16 of 25 random SOCPs at :math:`10^{-5}` instead of 25, in 12
times the time. Before 3.4.0, while refinement was on, every iteration was a
check, and each was followed by a conjugate gradient solve, to the solvers'
best accuracy, for the residual SCS computed, whether or not the solve had
missed its tolerance. Against that rule, in single precision at
:math:`10^{-6}` (100 random SOCPs, two starting points, two machines) 3
problems were lost and 6 gained, at about the same time (0.97 to 1.04 times
in shifted geometric mean); at :math:`10^{-5}` nothing changed but the time
(0.83 to 1.00 times). In double precision this refinement did not fire on the
problems measured.

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
  In single precision it tends to stay on once it turns on (see "Rounding
  level" above).

**Who computes the multiplies.** The multiplies by :math:`K` are made inside
the linear solve (:code:`scs_solve_lin_sys`, see :ref:`new_linear_solver`).
With the direct solvers SCS keeps each iteration's solution and passes it to
the next solve as its warm start; with the indirect solvers it asks for
refinement steps. Per iteration SCS asks for:

.. list-table::
   :header-rows: 1

   * - Iteration
     - Solves
     - Multiplies by :math:`K`
   * - Off, not a check
     - one plain solve
     - none
   * - Off, a check (every 25)
     - one plain solve, with its residual
     - 1
   * - On, not a check
     - one solve warm-started from the previous solution (correction form)
     - 1
   * - On, a check
     - the warm-started solve with its residual, and a plain solve with its
       residual
     - 3
   * - CG, off, a check
     - CG, with its residual
     - 1
   * - CG, on, not a check
     - CG with refinement steps: its residual and, where the solve missed its
       tolerance, up to two CG solves from zero for the residual
     - 1, and 1 per step
   * - CG, on, a check
     - the same, and a plain CG solve with its residual
     - 2, and 1 per step

By default the multiplies are SCS's own products with the problem data (the
multiply by :math:`A^\top` runs on OpenMP's threads when SCS is built with
OpenMP, and stays deterministic: each entry is summed by one thread in a fixed
order). MKL Pardiso computes them with MKL's threaded sparse matrix-vector
product instead, and cuDSS from 60,000 nonzeros in :math:`K` computes them on
the GPU, where the residual stays (the warm start is copied to the GPU). The
solve for :math:`(c, -b)` makes its multiplies the same way. The decisions
above are the same either way; only the cost of the multiply and its
rounding change.

What the solver's own multiplies buy is a cheaper iteration while the
correction form is on. With cuDSS on 10 large LPs, QPs and SOCPs at
:math:`10^{-8}` (NVIDIA GB10), such an iteration cost 1.19 times a plain one
with the multiplies on the GPU, against 1.48 with SCS's on the CPU (on BOYD1
1.61 ms against 2.63, a plain iteration taking 1.33; measured while the warm
start stayed on the GPU, see :ref:`cudss_solver` for the cost of copying
it). With MKL Pardiso on 33
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
stored nonzeros) gets SCS's own products instead; ``BLAS64=1`` (Meson's
``use_blas64``) provides MKL's at any size. If they cannot be set up
otherwise (out of memory), setup fails with a message. Built with
``make mkl MKL_NO_OWN_PRODUCTS=1`` (see :ref:`compile_flags`), MKL Pardiso
uses SCS's products at every size; that is slower (1.45 instead of 1.32
times a plain iteration while the correction form is on), and the CI uses it
so that this path stays tested.

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
:math:`2n` per vector while it runs (:math:`150n` by default). That is much
more than CG itself needs (a few vectors of length :math:`n`). Two constants
in :code:`include/glbopts.h` set it, and both are fixed when SCS is compiled:

* :code:`DEFLATE_VECTORS` (default 30) is the number of deflation vectors.
  More vectors take more of the smallest eigenvalues out of the solves, for
  more memory and more work per solve. :code:`0` turns deflation off, which
  gives plain preconditioned CG with no extra memory.
* :code:`EIGCG_WINDOW` (default 60) is the size of the eigCG window, the
  number of CG vectors that a harvest keeps in memory. A larger window
  compresses the CG history less often, for more memory while the harvest
  runs. The window is never smaller than three times the number of deflation
  vectors, so with the default 30 vectors it is 90 and a smaller
  :code:`EIGCG_WINDOW` has no effect unless :code:`DEFLATE_VECTORS` is lowered
  as well.

To change them, define them when compiling. With the Makefile, pass them in
:code:`CFLAGS` and rebuild from scratch (:code:`make purge` first, since
:code:`make` does not rebuild objects when only the flags change):

.. code:: bash

  make purge
  make CFLAGS="-DDEFLATE_VECTORS=0"

With CMake, pass them in :code:`CMAKE_C_FLAGS` when configuring:

.. code:: bash

  cmake -DCMAKE_C_FLAGS="-DDEFLATE_VECTORS=10 -DEIGCG_WINDOW=40" ..

A negative value stops the compilation with an error. Lowering
:code:`DEFLATE_VECTORS`, or setting it to 0, is the way to save memory with the
indirect solver on a problem where :math:`n` is very large. See also
:ref:`compile_flags`.

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
harvest eigenvectors or deflate, so :code:`DEFLATE_VECTORS` and
:code:`EIGCG_WINDOW` have no effect on it. Usually the direct GPU solver will be faster than this solver.

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
(``CUDSS_MIN_PRODUCT_NNZ`` in ``include/glbopts.h``; set it with
``make cudss CUDSS_MIN_PRODUCT_NNZ=...``, see :ref:`compile_flags`), the cuDSS
solver computes the multiplies by :math:`K` of the per-iteration refinement
(see :ref:`refinement`) on the GPU: one cuSPARSE sparse matrix-vector product
and one cuBLAS reduction. The residual stays on the GPU; the warm start, the
previous solution that SCS keeps on the host, is copied to the GPU once per
iteration with the correction form on, one vector of :math:`n + m` numbers
more than a plain solve copies. cuSPARSE has no symmetric product, so :math:`K`
is stored a second time with both triangles, about twice the nonzeros of the
copy cuDSS factorizes. A call costs about 0.04 ms whatever the size, so smaller
systems use SCS's product on the CPU, which was faster there. On an NVIDIA
GB10, per multiply, the GPU took 0.12 to 0.37 times the CPU's time above
:math:`3.7 \cdot 10^5` nonzeros (0.18 against 1.27 ms on BOYD1).
In the solves, on 10 large LPs, QPs and SOCPs at :math:`10^{-8}`, an
iteration with the correction form on cost 1.19 times a plain one, against
1.48 with the multiplies on the CPU, and the time per iteration fell by 4%.
That was measured with the warm start kept on the GPU. Copying it from the
host, as now, added 0.02 to 0.12 ms to such a solve on the same problems (2
to 9% of a plain solve): an estimated 0.14% of their total solve time, and at
most 3.3% on one (UBH1, with the correction form on for 64% of its
iterations). On a discrete GPU, where the copy crosses PCIe, it has not been
measured.
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

A linear solver implements the struct :code:`ScsLinSysWork`, which holds its
workspace, and the five functions in :code:`include/linsys.h` (documented
below); the :code:`linsys` directory has seven worked examples. Four are
simple: set up (for instance factorize :math:`K`), update for a new
:math:`R`, free, and name. The fifth, :code:`scs_solve_lin_sys`, solves
:math:`K x = b` and also takes a warm start, a tolerance and a number of
refinement steps, and reports the residual :math:`\|b - K x\|_\infty`:
SCS's :ref:`iterative_refinement` runs through it. The solvers in
:code:`linsys` write only their plain solve and get the rest from a shared
routine (:code:`linsys/scs_matrix.h`). A new solver built with SCS's sources
does the same in five steps:

1. :code:`#include "scs_matrix.h"` and put a :code:`ScsLinSysShared sh;`
   member in your :code:`ScsLinSysWork`.

2. Write your plain solve, which overwrites :code:`b` with :math:`K^{-1} b`
   and returns 0 (nonzero on failure). A direct solver ignores :code:`s` and
   :code:`tol`; an iterative one starts from :code:`s` (its first :math:`n`
   entries; from zero when it is :code:`SCS_NULL`) and stops at
   :code:`tol`::

       static scs_int plain_solve(ScsLinSysWork *w, scs_float *b,
                                  const scs_float *s, scs_float tol);

3. At the end of :code:`scs_init_lin_sys_work`, after your own setup::

       if (SCS(lin_sys_shared_init)(&w->sh, A, P, diag_r, 0, plain_solve,
                                    SCS_NULL, SCS_NULL) != 0) {
         scs_free_lin_sys_work(w);
         return SCS_NULL;
       }

   with 1 in place of the first 0 for an iterative solver such as conjugate
   gradients (its warm start is then the solver's starting point, and each
   refinement step a solve from zero of the residual, computed in the full
   space; see :code:`linsys.h`). This allocates three vectors of
   :math:`n + m` numbers and keeps the three pointers: SCS keeps :code:`A`,
   :code:`P` and :code:`diag_r` alive until :code:`scs_free_lin_sys_work`,
   leaves :code:`A` and :code:`P` unchanged and writes each new :math:`R`
   into the same :code:`diag_r`, and the residuals are computed from them.

4. Implement :code:`scs_solve_lin_sys` as one call::

       scs_int scs_solve_lin_sys(ScsLinSysWork *w, scs_float *b,
                                 const scs_float *s, scs_float tol,
                                 scs_int max_refine_steps, scs_float *res) {
         return SCS(lin_sys_shared_solve)(&w->sh, w, b, s, tol,
                                          max_refine_steps, res);
       }

5. Call :code:`SCS(lin_sys_shared_free)(&w->sh)` in
   :code:`scs_free_lin_sys_work` (it is safe on a zeroed :code:`sh`).

With that, SCS's refinement behaves with your solver as with the bundled
ones: the warm start (the correction form and its comparison with zero, or
the starting point of an iterative solver), the refinement steps, and the
residuals are the shared routine's, with SCS's own
multiplies by :math:`K` from the problem data. Only your plain solve's
accuracy and cost differ. :code:`linsys/scs_matrix.h` is not installed with
SCS's headers, so this is for solvers compiled with SCS's sources, as every
solver outside this repository we know of is.

**Your own multiplies (optional).** A solver that holds :math:`K` anyway can
compute the multiplies with its own copy instead, for instance with a
threaded sparse matrix-vector product, by passing a residual function in
place of the first :code:`SCS_NULL` in step 3::

    static scs_int my_residual(ScsLinSysWork *w, const scs_float *b,
                               const scs_float *x, scs_float *r,
                               scs_float *nr, scs_float *nb);

It computes :math:`r = b - K x` with your :math:`K` (as of the last update of
:math:`R`) into :code:`r`, :math:`\|r\|_\infty` in :code:`*nr` and
:math:`\|b\|_\infty` in :code:`*nb` (unless :code:`nb` is
:code:`SCS_NULL`), and returns 0 (nonzero on failure, which fails the solve
with a message asking to open an issue). MKL Pardiso does this with MKL's
sparse BLAS (:code:`prod_residual` in :code:`linsys/mkl/direct/private.c`). A
solver whose vectors live elsewhere can also pass a correct function in place
of the second :code:`SCS_NULL`::

    static scs_int my_correct(ScsLinSysWork *w, const scs_float *x,
                              const scs_float *r, scs_float *out);

It sets :code:`out` to :math:`x + K^{-1} r`, where :code:`r` is that of the
residual call the shared routine always makes just before and :code:`x` has
the values of that call's :code:`x`. So the residual function may keep
:math:`r` where the solver keeps its vectors, and need not write the
:code:`r` it is given; if it keeps :math:`x` too, it keeps a copy, not the
pointer: :code:`x` may then be another vector with the same values, and
:code:`out` the very vector the residual call was given as :code:`x`. cuDSS does
this on the GPU (:code:`dev_residual` and :code:`dev_correct` in
:code:`linsys/cudss/direct/private.c`), so that the residual never leaves
the device. Decide whether to use your own multiplies once, in
:code:`scs_init_lin_sys_work`, by a rule you state up front (MKL: always,
except a :math:`K` too large for 32-bit MKL integers in a ``DLONG`` build
without ``BLAS64``; cuDSS: from :code:`CUDSS_MIN_PRODUCT_NNZ` nonzeros),
never by whether memory or a library happened to be available: if they
should be used but cannot be set up, make :code:`scs_init_lin_sys_work` fail
with a message saying why.

**Without the shared routine.** A solver may implement
:code:`scs_solve_lin_sys` entirely itself, but then it must honour the whole
of :code:`linsys.h`. :code:`res` must be the true residual of the returned
:math:`x`, because SCS turns its refinement on and off by comparing it with
its own residuals. A direct solver given :code:`s` must return the
correction form from it, and solve plainly when :math:`\|b - K s\|_\infty
\ge \|b\|_\infty`. With :code:`max_refine_steps` above 0 it must take at most
that many steps :math:`x + K^{-1}(b - K x)` (a CG solver too, with the
residual computed in the full space and solved from zero), keep a step only
if it lowers the residual, and stop after one that does not halve it. SCS
asks for such steps to solve for :math:`(c, -b)` (see :ref:`g_refinement`),
and stops with an error if that solve fails. A direct solver that ignores
:code:`s` still gives correct answers, but SCS's refinement then cannot make
its solves more accurate, which on badly scaled problems can keep SCS from
reaching tight tolerances.

**Tests.** :code:`test_lin_sys_solve` in the test suite checks the contract
against residuals it computes itself: warm starts near and far from the
solution (the correction form and its fallback), refinement and its limits
(steps kept only when better, the halving stop, :code:`max_refine_steps`),
and the residual with the :math:`K` of the solver's rule,
before and after an update of :math:`R`. :code:`test_refine_perturbed` makes
the solver factorize a different :math:`R` from the one SCS measures with, so
that SCS's refinement fires (and, with a mismatch too large to contract, its
safeguard); a solver that measures with its own :math:`K` follows the update,
sees no mismatch, and is skipped there with a message (MKL Pardiso, unless
built with ``MKL_NO_OWN_PRODUCTS``; cuDSS from :code:`CUDSS_MIN_PRODUCT_NNZ`
nonzeros). :code:`test_g_accuracy` checks that the solve for
:math:`(c, -b)` reaches machine precision on a badly scaled problem, and
that a solver given a far different :math:`R` (skipped likewise) makes SCS
stop with :code:`SCS_FAILED` and the message asking to open an issue.

**From SCS 3.3.** The solve was :code:`scs_solve_lin_sys(w, b, s, tol)`.
:code:`max_refine_steps` and :code:`res` are new, so a solver written for 3.3
no longer compiles. Rename its solve :code:`plain_solve` (it keeps its 3.3
meaning) and follow the steps above; the solver then behaves as the bundled
solvers do. A solver written against the unreleased development versions
after 3.3, which had :code:`scs_lin_sys_residual`, :code:`scs_lin_sys_keep` and
:code:`scs_lin_sys_correct`, deletes them; its own multiplies, if it has
them, become a residual function as above.

.. doxygentypedef:: ScsLinSysWork

.. doxygenfile:: linsys.h

