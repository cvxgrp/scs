.. This chapter is a first draft written for maintainer review. Every
   recommendation below should be checked by someone who knows the solver
   better than its documentation does before it ships.

.. _best_practices:

Best practices
==============

SCS is a first-order method. That gives it its strengths (it scales, it
detects infeasibility, it warm-starts, and it handles every cone the same way)
and its one weakness: it converges quickly to modest accuracy and slowly to
high accuracy. Most of the advice below follows from that single fact.

Scale your problem
------------------

The single most effective thing you can do for SCS is to hand it well-scaled
data. Aim for the nonzero entries of :math:`A`, :math:`P`, :math:`b` and
:math:`c` to be within a few orders of magnitude of one another, ideally around
one. Choose units for your variables and constraints with that in mind: a
constraint expressed in millimetres next to one expressed in kilometres will
hurt.

SCS :ref:`equilibrates <equilibration>` the data before solving (the
:code:`normalize` setting, on by default) and adapts its internal metric during
the solve (:code:`adaptive_scale` and :code:`adaptive_diag_scale`, both on by
default). These fix moderate scaling problems automatically, and you should
leave them on. They cannot fix extreme ones: if your data spans ten orders of
magnitude, rescale it yourself first.

If you have turned :code:`normalize` off for some reason, the :code:`scale`
setting becomes the main knob; values between :code:`0.01` and :code:`10`
cover most problems.

Choose tolerances deliberately
------------------------------

The stopping criteria are relative (see :ref:`termination <termination>`), so
:code:`eps_abs` and :code:`eps_rel` are dimensionless. The defaults of
:code:`1e-4` give a solution that is accurate to roughly four digits in the
residuals, which is enough for most applications, and SCS typically reaches it
in a few hundred to a few thousand iterations.

Tighter tolerances cost more than you might expect. Each additional decade of
accuracy can cost a multiple of the iterations of the last one, particularly on
degenerate problems where the optimal set is not a single point. Before
tightening, ask whether you need the extra digits in the solution or only a
smaller objective error: the objective converges faster than the iterates. If
you do need :code:`1e-8` or better on a small problem, an interior-point solver
may be the better tool; see :ref:`when_not_scs` below.

:code:`eps_infeas` (default :code:`1e-7`) controls how readily SCS declares a
problem infeasible or unbounded. Decrease it if you see spurious certificates on
a problem you know to be feasible; see :ref:`infeasibility`.

Solve sequences of problems with a cached workspace
---------------------------------------------------

If you solve many problems that share :math:`A` and :math:`P` and differ only
in :math:`b` and :math:`c`, initialize once and reuse the workspace. In C that
is :code:`scs_init` followed by repeated :code:`scs_update` and
:code:`scs_solve` calls (see the :ref:`C interface <c_interface>`); in Python it
is one :code:`scs.SCS` object with :code:`update` and :code:`solve`. The matrix
factorization and the equilibration are computed once and reused, which is
usually the dominant cost of a solve.

Warm-start whenever the previous solution is close to the next one, as in a
model-predictive-control loop or when tracing a regularization path. Pass the
previous :math:`x`, :math:`y` and :math:`s`; SCS starts from them and
overwrites them with the new solution. See :ref:`warm_start` and the
:ref:`MPC example <py_mpc>`.

Pick the right linear solver
----------------------------

Every iteration solves one linear system with a fixed matrix, and the choice
of how to solve it is the main performance decision. The
:ref:`linear solver <linear_solver>` page describes each backend; as a rule of
thumb:

* **Default (sparse direct, QDLDL).** Factorizes once and reuses the factors.
  Right for most problems up to a few hundred thousand nonzeros. Always
  available.
* **MKL Pardiso.** A parallel direct solver, often several times faster than
  the default on large problems. Bundled in the x86-64 Linux and Windows Python
  wheels and selected automatically there; see :ref:`mkl`.
* **Apple Accelerate.** Available on macOS; selected explicitly. See
  :ref:`apple_accelerate`.
* **Sparse indirect (conjugate gradient).** Never factorizes, so it handles
  problems too large to factorize and can be run matrix-free. Each iteration
  is cheaper but there are more of them, and it is less predictable. See
  :ref:`indirect`.
* **cuDSS on the GPU.** The recommended GPU backend, a direct solver on the
  device. The older GPU indirect solver is a legacy backend that only pays off
  for very large problems. See :ref:`cudss_solver` and :ref:`gpu_indirect`.

For the Python interface, :code:`linear_solver=scs.LinearSolver.AUTO` (the
default) picks the best available direct solver for the platform.

Leave acceleration on, but know how to turn it off
--------------------------------------------------

:ref:`Anderson acceleration <acceleration>` is on by default (type-I, memory
:code:`10`, applied every :code:`5` iterations) and usually reduces the
iteration count substantially, especially at tight tolerances. It can also
destabilize the solve on some problems. If you see the residuals jump around
rather than decrease, or a solve that stalls, set
:code:`acceleration_lookback=0` and compare. The :code:`aa_stats` field of the
returned info reports how many accelerated steps were accepted and why any were
rejected, which is the place to look when tuning
:code:`acceleration_regularization`.

Read the log
------------

With :code:`verbose` on, SCS prints the problem dimensions, the cones, the
settings and the linear solver, then one line every few iterations with the
primal residual, dual residual, duality gap, objective, current :code:`scale`
and elapsed time, and finally a status line and the timings broken down into
linear-system solves, cone projections and acceleration. Two habits pay off:

* Check the header the first time you solve a new problem. If the dimensions
  or cones are not what you expect, the data is wrong and nothing downstream
  will help.
* Watch which residual is lagging. If one of primal and dual residual is far
  ahead of the other, the problem is badly scaled between constraints and
  objective; the adaptive scaling will try to correct it, and you can help by
  rescaling.

The status :code:`solved_inaccurate` means SCS stopped at :code:`max_iters` or
:code:`time_limit_secs` before meeting the tolerances; the returned point is
the best found and is often perfectly usable. See the :doc:`help page </help/index>` for what to
try next.

Interpret infeasibility certificates correctly
----------------------------------------------

SCS returns a certificate rather than a solution when it finds one (see
:ref:`infeasibility`). Before trusting a certificate you did not expect, check
the two most common causes of a spurious one: the rows of :math:`A` and
:math:`b` not following the required :ref:`cone order <cones>`, and a sign or
transposition error in the data. Setting :code:`write_data_filename` dumps
exactly what SCS received, which is the fastest way to inspect it.

Differentiating through SCS
---------------------------

Packages such as diffcp and cvxpylayers compute derivatives of the solution
with respect to the problem data using the solution SCS returns. The accuracy
of those derivatives is limited by the accuracy of that solution, and default
tolerances are usually too loose for the purpose: pass tight tolerances (of
order :code:`1e-8` or better) when a solve feeds a derivative computation, and
expect finite-difference checks against such derivatives to be limited by the
solve accuracy divided by the step size.

.. _when_not_scs:

When to use something else
--------------------------

SCS is the right tool for large problems, for problems with semidefinite,
exponential or power cones, for problems that need infeasibility detection,
and for sequences of related problems. It is not the best tool for a small
problem that needs a highly accurate answer: an interior-point solver will
reach :code:`1e-8` in tens of iterations where a first-order method needs
thousands. For pure quadratic programs at modest accuracy a dedicated QP solver
may also be faster. Modeling layers such as CVXPY choose a default solver by
problem class for exactly this reason, and SCS is their default for the classes
where it is strongest.

Reporting problems
------------------

If SCS misbehaves on a problem, set :code:`write_data_filename`, attach the
resulting file to a GitHub issue together with the settings you used and the
solver log, and say what you expected instead. That is almost always enough to
reproduce the behavior; see the :doc:`help page </help/index>`.
