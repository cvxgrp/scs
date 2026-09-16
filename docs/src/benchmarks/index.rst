.. _benchmarks:

Benchmarks
==========

This page compares SCS 3.3 against other open-source solvers on standard
quadratic program (QP), linear program (LP) and semidefinite program (SDP)
test sets. All results were produced with the open
`solver_benchmarks <https://github.com/bodono/solver_benchmarks>`_ harness,
which feeds every solver the same problem in its native form, applies the same
time limit, and then checks every returned solution independently. The raw
results, the campaign configuration and the plotting scripts live in that
repository, so every number here can be regenerated.

The short version: on the QP and LP test sets SCS is competitive with the best
interior-point codes and clearly ahead of the other first-order solvers, the
cuDSS GPU backend is the fastest solver we tested on the largest quarter of the
QP problems, and on the large LPs of the Mittelmann set SCS with cuDSS solves
more instances than any other solver. On the SDP test sets the interior-point solvers are faster on
the small, ill-conditioned control and truss instances, while SCS is the
fastest solver on the large sparse combinatorial relaxations; see
:ref:`bench_sdp` before choosing a solver for semidefinite problems.

.. _bench_method:

Methodology
-----------

**Verification.** A solver's own "optimal" status is not taken at face value.
For every returned solution the harness recomputes the relative primal
residual, dual residual and duality gap from the problem data, and a solve
counts as a success only when all three are below ten times the requested
tolerance (so :math:`10^{-3}` for the :math:`10^{-4}` runs and :math:`10^{-5}`
for the :math:`10^{-6}` runs). A solve that a solver reports as optimal but
whose residuals fail this check is counted as a failure; a solve reported as
inaccurate or timed out whose residuals pass is counted as a success. This
matters because solvers define "relative tolerance" differently, and because
first-order methods, SCS included, occasionally declare convergence on badly
scaled problems where the unscaled residuals are large. Applying the same
independent check to every solver makes the comparison tolerance-fair.

**Tolerances.** The first-order solvers (SCS, OSQP, ProxQP, PDLP and cuOpt)
were run at a requested tolerance of :math:`10^{-4}` for the headline plots
and :math:`10^{-6}` for the high-accuracy plots, passing the value through
each solver's own absolute and relative settings; for SCS this means
``eps_abs = eps_rel = tol`` with the iteration limit raised so that only the
time limit can stop it. For reference, the solvers' own defaults are
:math:`10^{-4}` for SCS and cuOpt, :math:`10^{-3}` for OSQP, :math:`10^{-5}`
(absolute only) for ProxQP and :math:`10^{-6}` for OR-Tools PDLP, each in its
own measure of the residual. The
interior-point and simplex solvers (Clarabel, PIQP, HiGHS, SDPA and CVXOPT)
are shown from their :math:`10^{-6}` runs in every plot. Nobody runs an interior-point solver at :math:`10^{-4}`: their defaults
are :math:`10^{-7}` to :math:`10^{-8}`, the last few iterations are nearly
free, and a loose setting would show them faster than a user ever sees them.
It also removes a real comparability problem: Clarabel's termination test can
leave large unscaled residuals when some terms in it are very large, and at
:math:`10^{-4}` on the SDP sets its solutions were up to a hundred times
looser than everyone else's. Every solve, whatever was requested, is then
verified independently at ten times the plot's tolerance. The first-order
solvers were also run at :math:`10^{-5}`; see :ref:`bench_sens` for how
much that changes the picture.

A requested tolerance means different things to different solvers, so the
table below reports the accuracy actually achieved, measured the same way for
everyone: the largest of the three relative KKT residuals over the solves each
solver itself reported optimal, at the :math:`10^{-4}` setting. The
first-order solvers land at the requested tolerance; the interior-point and
simplex codes overshoot it by orders of magnitude, except Clarabel on SDP;
cuOpt's LP path and ProxQP have a long tail of returns that are far from
optimal.

.. list-table:: Achieved accuracy: median and 90th percentile of the largest relative KKT residual over the solves each solver reported optimal, at the 1e-4 and 1e-6 settings
   :header-rows: 2
   :widths: 22 13 13 13 13 13 13

   * - Solver
     - QP 1e-4
     -
     - LP 1e-4
     -
     - SDP 1e-4
     -
   * -
     - median
     - 90th
     - median
     - 90th
     - median
     - 90th
   * - SCS (CPU)
     - 1e-5
     - 1e-4
     - 3e-5
     - 9e-5
     - 9e-5
     - 1e-4
   * - SCS (GPU, cuDSS)
     - 2e-5
     - 9e-5
     - 3e-5
     - 1e-4
     - 8e-5
     - 1e-4
   * - OSQP
     - 3e-5
     - 1e-4
     - 6e-5
     - 1e-4
     - --
     - --
   * - PDLP (OR-Tools)
     - --
     - --
     - 1e-4
     - 8e-4
     - --
     - --
   * - ProxQP
     - 8e-5
     - 3e-2
     - --
     - --
     - --
     - --
   * - cuOpt (GPU)
     - 2e-8
     - 1e-6
     - 1e-4
     - 1e+00
     - --
     - --
   * - Clarabel
     - 1e-5
     - 6e-5
     - 2e-5
     - 1e-4
     - 3e-4
     - 5e-3
   * - PIQP
     - 4e-10
     - 1e-7
     - 2e-10
     - 6e-8
     - --
     - --
   * - HiGHS
     - 1e-7
     - 3e-4
     - 8e-16
     - 1e-10
     - --
     - --
   * - CVXOPT
     - --
     - --
     - --
     - --
     - 2e-5
     - 8e-5

**Timing.** The reported time is wall-clock time for the solver call, including
any presolve, factorization and GPU transfer, but excluding reading the problem
from disk. QP and LP solves were limited to 300 s and SDP solves to 900 s, passed to each
solver through its own time-limit setting; since wall time also includes
setup and reading, a solve counts only if it finished within the limit plus a
60 s grace, the same for every solver, and the harness killed the worker at
that point regardless of what the solver reported. A
performance profile shows, for each solver, the fraction of problems solved
within a factor :math:`\tau` of the fastest solver on that problem; failures
never count as solved. The shifted geometric mean uses a shift of 10 s and
charges each failure 1000 s. Times below 10 ms are floored at 10 ms before
computing ratios.

**Hardware.** All CPU solvers ran in identical 4-core Linux x86-64 containers
(Modal), on which the ``scs`` wheel selects the MKL Pardiso linear system
solver by default. The two GPU solvers, SCS with cuDSS and NVIDIA cuOpt, ran on
an NVIDIA A100 80GB with 8 host cores. GPU results include host-device transfer
and cuDSS analysis time, so small problems pay a fixed overhead of roughly half
a second.

**Solvers.** SCS 3.3.1 (CPU with MKL Pardiso, and GPU with cuDSS),
Clarabel, PIQP, OSQP, ProxQP, HiGHS, PDLP (OR-Tools 9.15), NVIDIA cuOpt 26.8,
CVXOPT and SDPA, each at the latest release on PyPI at the time of the run.
Commercial solvers were not included.

**Problem sets.** QP: Maros-Meszaros (138) and the QPLIB continuous convex
subset (19). LP: Netlib (feasible), Kennington and the MIPLIB 2017
LP relaxations up to 20 MB; the Mittelmann LP set is its own section. SDP: SDPLIB and the
Mittelmann SDP set. "Largest quartile" means the quarter of each family with
the most nonzeros in the constraint matrix (plus the Hessian for QPs).

**Exclusions.** Clarabel ran out of memory (64 GB) on the SDPLIB instances
with PSD blocks of order 500 or more; those instances are counted as failures
for it. SDPA does not take a time limit but never needed one. The SDPLIB archive's
``maxG55`` and ``maxG60`` files are corrupt and were dropped for all solvers.
cuOpt's QP path is an interior-point method whose factorization failed with a
numerical error on a subset of the Maros-Meszaros problems; those count as
failures.

.. _bench_qp:

Quadratic programs
------------------

157 problems: Maros-Meszaros (138) and the convex continuous QPLIB instances
(19). On the largest quarter of these, SCS with cuDSS
has the lowest shifted geometric mean solve time of any solver, SCS on the CPU
sits between Clarabel and PIQP, and the two SCS variants verify the most
solutions. Over all 157 problems the interior-point solvers are faster on the
small instances, where an SCS solve is dominated by fixed setup cost, but SCS
solves nearly as many problems as they do.

.. figure:: ../files/bench/qp_1e-4_profile_largest.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/qp_1e-4_geomean_largest.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/qp_1e-4_profile.png
   :width: 90 %
   :align: center

At the tighter :math:`10^{-6}` tolerance the interior-point solvers pull ahead,
as expected for a first-order method, but SCS still verifies more solutions
than every solver other than PIQP and Clarabel, and remains the fastest
first-order solver by a wide margin.

.. figure:: ../files/bench/qp_1e-6_profile_largest.png
   :width: 90 %
   :align: center

.. list-table:: QP: verified solves and shifted geometric mean time (s); all 157 problems / largest quartile (40)
   :header-rows: 1
   :widths: 26 12 12 12 12 12 12 12 12

   * - Solver
     - solved 1e-4
     - gm 1e-4
     - solved 1e-4 (largest)
     - gm 1e-4 (largest)
     - solved 1e-6
     - gm 1e-6
     - solved 1e-6 (largest)
     - gm 1e-6 (largest)
   * - Clarabel
     - 152
     - 3.2
     - 37
     - 12.8
     - 148
     - 4.5
     - 34
     - 20.0
   * - PIQP
     - 150
     - 3.7
     - 34
     - 18.7
     - 150
     - 3.7
     - 34
     - 18.7
   * - SCS (GPU, cuDSS)
     - 146
     - 6.4
     - 38
     - 9.5
     - 140
     - 12.5
     - 33
     - 26.8
   * - SCS (CPU, MKL Pardiso)
     - 144
     - 7.6
     - 38
     - 12.9
     - 138
     - 11.7
     - 33
     - 29.6
   * - OSQP
     - 138
     - 11.1
     - 32
     - 30.2
     - 121
     - 28.1
     - 30
     - 43.2
   * - cuOpt (GPU)
     - 109
     - 35.2
     - 32
     - 23.1
     - 103
     - 41.4
     - 27
     - 38.7
   * - HiGHS
     - 100
     - 50.2
     - 13
     - 288.3
     - 84
     - 83.5
     - 11
     - 349.8
   * - ProxQP
     - 81
     - 106.3
     - 10
     - 439.7
     - 84
     - 93.8
     - 12
     - 375.1

.. _bench_lp:

Linear programs
---------------

344 problems: Netlib, Kennington and the MIPLIB 2017 LP relaxations up to
20 MB. SCS is not an LP solver and is not marketed as one,
so this comparison is included mainly to show that it is not a bad one: on the
largest quarter of the LP set SCS with cuDSS verifies as many solutions as
HiGHS (dual simplex) and is second only to it in geometric mean time, ahead of
the interior-point solvers PIQP and Clarabel. Both PDLP implementations, which
are first-order LP methods, trail SCS by a wide margin under independent
verification (see the notes below).

.. figure:: ../files/bench/lp_1e-4_profile_largest.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/lp_1e-4_geomean_largest.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/lp_1e-4_profile.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/lp_1e-6_profile_largest.png
   :width: 90 %
   :align: center

.. list-table:: LP: verified solves and shifted geometric mean time (s); all 344 problems / largest quartile (86)
   :header-rows: 1
   :widths: 26 12 12 12 12 12 12 12 12

   * - Solver
     - solved 1e-4
     - gm 1e-4
     - solved 1e-4 (largest)
     - gm 1e-4 (largest)
     - solved 1e-6
     - gm 1e-6
     - solved 1e-6 (largest)
     - gm 1e-6 (largest)
   * - HiGHS
     - 331
     - 4.6
     - 74
     - 25.5
     - 331
     - 4.6
     - 74
     - 25.5
   * - PIQP
     - 315
     - 8.5
     - 65
     - 43.2
     - 309
     - 9.7
     - 64
     - 44.0
   * - Clarabel
     - 323
     - 9.5
     - 70
     - 60.7
     - 311
     - 12.4
     - 65
     - 74.6
   * - SCS (GPU, cuDSS)
     - 315
     - 10.8
     - 74
     - 34.6
     - 298
     - 16.5
     - 66
     - 50.1
   * - SCS (CPU, MKL Pardiso)
     - 314
     - 10.9
     - 71
     - 49.4
     - 294
     - 18.8
     - 61
     - 81.3
   * - PDLP (OR-Tools)
     - 280
     - 24.0
     - 47
     - 144.8
     - 172
     - 130.7
     - 39
     - 234.6
   * - OSQP
     - 256
     - 41.0
     - 44
     - 210.2
     - 176
     - 131.9
     - 23
     - 480.3
   * - cuOpt (GPU)
     - 133
     - 163.7
     - 24
     - 280.1
     - 123
     - 187.4
     - 20
     - 342.0

.. _bench_sdp:

Semidefinite programs
---------------------

98 problems: SDPLIB and the Mittelmann SDP set. This is the family where an
interior-point method is usually the right choice, and the plots say so: SDPA
and CVXOPT solve the small, ill-conditioned ``control``, ``truss``, ``arch``
and ``gpp`` instances in seconds where SCS needs hundreds of thousands of
iterations and often hits the 900 s limit. SCS is the fastest solver on the
large sparse combinatorial relaxations (``theta``, ``mcp``, ``maxG``, ``qpG``
and ``equalG``), where the interior-point methods either run out of time or,
in Clarabel's case, cannot form the dense scaling block at all. The GPU does
not help on SDPs: the time is spent in the eigendecompositions of the cone
projection, not in the linear system.

If your problem has a few large PSD blocks and moderate accuracy is enough,
SCS is a good choice; if it has many small blocks or is badly conditioned, use
an interior-point solver.

.. figure:: ../files/bench/sdp_1e-4_profile.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/sdp_1e-4_geomean.png
   :width: 90 %
   :align: center

.. list-table:: SDP: verified solves and shifted geometric mean time (s); all 98 problems / largest quartile (28)
   :header-rows: 1
   :widths: 26 12 12 12 12 12 12 12 12

   * - Solver
     - solved 1e-4
     - gm 1e-4
     - solved 1e-4 (largest)
     - gm 1e-4 (largest)
     - solved 1e-6
     - gm 1e-6
     - solved 1e-6 (largest)
     - gm 1e-6 (largest)
   * - SDPA
     - 72
     - 38.0
     - 23
     - 35.5
     - 66
     - 53.5
     - 23
     - 35.5
   * - CVXOPT
     - 73
     - 51.6
     - 22
     - 62.4
     - 59
     - 108.0
     - 20
     - 88.4
   * - Clarabel
     - 61
     - 80.5
     - 16
     - 140.3
     - 47
     - 139.1
     - 11
     - 240.2
   * - SCS (CPU, MKL Pardiso)
     - 76
     - 98.0
     - 18
     - 674.2
     - 50
     - 192.1
     - 2
     - 989.5
   * - SCS (GPU, cuDSS)
     - 72
     - 126.6
     - 13
     - 832.2
     - --
     - --
     - --
     - --

.. _bench_lpbig:

Large linear programs: the Mittelmann set
-----------------------------------------

The LP test sets above are dominated by small and medium instances, so we also
ran the 37 problems of `Hans Mittelmann's LP benchmark set
<https://plato.asu.edu/ftp/lptestset/>`_, the standard collection of large,
hard LPs: between 100,000 and 126 million nonzeros, with several instances of
10 to 40 million variables. Every solver ran at tolerance :math:`10^{-4}`
with an 1800 s limit in 64 GB containers (4 cores, or an A100 80GB for the two
GPU solvers). Mittelmann's own runs allow several hours per instance and use
faster machines, so the simplex and interior-point codes time out here far
more often than they do in his tables; the point of this set for us is the
size of the problems, not a re-run of his benchmark.

This is where the cuDSS backend pays off. SCS on the GPU verifies more
solutions than any other solver and has by far the lowest geometric mean
time, and SCS on the CPU is second. The two other first-order codes, PDLP and
cuOpt, are the natural comparison: OR-Tools PDLP verifies about two thirds as many
solutions as SCS with cuDSS, and cuOpt's PDLP, although it reports almost
every instance optimal, mostly fails the independent residual check (see the
notes below).

.. figure:: ../files/bench/lpbig_1e-4_profile.png
   :width: 90 %
   :align: center

.. figure:: ../files/bench/lpbig_1e-4_geomean.png
   :width: 90 %
   :align: center

.. list-table:: Mittelmann LP set: verified solves out of 37 and shifted geometric mean time (s), tolerance 1e-4, 1800 s limit
   :header-rows: 1
   :widths: 40 20 20

   * - Solver
     - verified solves
     - geometric mean (s)
   * - SCS (GPU, cuDSS)
     - 23
     - 202
   * - SCS (CPU, MKL Pardiso)
     - 20
     - 316
   * - Clarabel
     - 20
     - 372
   * - PDLP (OR-Tools)
     - 14
     - 425
   * - PIQP
     - 12
     - 428
   * - cuOpt (GPU)
     - 4
     - 675
   * - HiGHS
     - 9
     - 804

Set-specific exclusions: Clarabel could not attempt ``L1_sixm250obs`` and
``L1_sixm1000obs`` within 64 GB (counted as failures); OR-Tools PDLP cannot
load ``Dual2_5000`` and ``dlr2`` because the model exceeds the 2 GB protobuf
limit (counted as failures); SCS with cuDSS ran out of GPU memory on
``thk_48`` (counted as a failure). Four instances (``bdry2``, ``Linf_520c``
and the two ``L1_sixm`` problems) are distributed in Netlib's compressed EMPS
format and were decoded with ``emps`` before use.

.. _bench_sens:

Sensitivity to the requested tolerance
--------------------------------------

The headline plots ask the first-order solvers for :math:`10^{-4}`. To show
what a tighter request costs, every first-order solver was also run at
:math:`10^{-5}` on the QP, LP and Mittelmann sets, verified at
:math:`10^{-4}`; the interior-point rows are unchanged, since their runs are
already tighter than that. The three landing-page rows at :math:`10^{-5}`:

.. figure:: ../files/bench/landing_grid_1e-5.png
   :width: 100 %
   :align: center

The LP picture barely moves and the Mittelmann set is still led by SCS with
cuDSS and SCS on the CPU, with a smaller margin. The QP row is where the
tolerance matters: SCS loses a few of the largest instances to the time limit
and its geometric mean roughly doubles, so at :math:`10^{-5}` Clarabel and
PIQP are faster on the largest QPs while SCS still matches Clarabel on the
number solved. The other first-order solvers lose more than SCS does from the
tighter request, OSQP and PDLP a quarter to a third of their large solves, so
SCS's margin over them widens. cuOpt's LP returns pass verification far more
often at :math:`10^{-5}` (219 of 344 against 133), though it remains slow.

.. list-table:: Verified solves and shifted geometric mean time (s) at the 1e-4 and 1e-5 settings; interior-point solvers unchanged (shown from their tightest run in both)
   :header-rows: 1
   :widths: 30 24 12 12 12 12

   * - Problem set
     - Solver
     - solved 1e-4
     - time 1e-4
     - solved 1e-5
     - time 1e-5
   * - Maros-Meszaros and QPLIB QPs, largest quartile (40)
     - Clarabel
     - 37
     - 13
     - 37
     - 13
   * - 
     - PIQP
     - 34
     - 19
     - 34
     - 19
   * - 
     - SCS (GPU, cuDSS)
     - 38
     - 9
     - 33
     - 21
   * - 
     - SCS (CPU, MKL Pardiso)
     - 38
     - 13
     - 35
     - 21
   * - 
     - cuOpt (GPU)
     - 32
     - 23
     - 31
     - 24
   * - 
     - OSQP
     - 32
     - 30
     - 31
     - 39
   * - 
     - HiGHS
     - 13
     - 288
     - 12
     - 313
   * - 
     - ProxQP
     - 10
     - 440
     - 12
     - 373
   * - Kennington and MIPLIB-relaxation LPs, largest quartile (86)
     - HiGHS
     - 74
     - 25
     - 74
     - 25
   * - 
     - SCS (GPU, cuDSS)
     - 74
     - 35
     - 72
     - 37
   * - 
     - PIQP
     - 65
     - 43
     - 65
     - 43
   * - 
     - SCS (CPU, MKL Pardiso)
     - 71
     - 49
     - 69
     - 57
   * - 
     - Clarabel
     - 70
     - 61
     - 66
     - 72
   * - 
     - PDLP (OR-Tools)
     - 47
     - 145
     - 41
     - 201
   * - 
     - cuOpt (GPU)
     - 24
     - 280
     - 24
     - 279
   * - 
     - OSQP
     - 44
     - 210
     - 34
     - 360
   * - Mittelmann LP set (37)
     - SCS (GPU, cuDSS)
     - 28
     - 135
     - 23
     - 202
   * - 
     - SCS (CPU, MKL Pardiso)
     - 24
     - 235
     - 20
     - 316
   * - 
     - Clarabel
     - 20
     - 372
     - 20
     - 372
   * - 
     - PDLP (OR-Tools)
     - 19
     - 325
     - 14
     - 425
   * - 
     - PIQP
     - 12
     - 428
     - 12
     - 428
   * - 
     - cuOpt (GPU)
     - 7
     - 528
     - 4
     - 675
   * - 
     - HiGHS
     - 9
     - 804
     - 9
     - 804

Notes on individual solvers
---------------------------

* **SCS.** A few SCS "solved" returns fail the independent residual check
  and count as failures above. Two of them are real false certificates
  (``QGROW7`` and ``QGROW22`` from Maros-Meszaros, and the MIPLIB relaxation
  ``neos-4413714-turia``): on these badly scaled problems the diagonal
  rescaling drives the internal scale to its floor, and the scaled residuals
  SCS monitors no longer track the unscaled ones.
* **cuOpt.** cuOpt's LP path is PDLP with an L2-relative stopping rule, and on
  problems with many variable bounds (Kennington, MIPLIB) its returned duals
  often have large infinity-norm stationarity residuals even though cuOpt
  reports the solve as optimal, and even though its own reported absolute
  dual residual is between :math:`3\times10^{1}` and :math:`3\times10^{2}`. Under the uniform check used here
  those solves count as failures; by its own status cuOpt reports 233 of the
  344 LPs optimal at :math:`10^{-4}`, and 36 of the 37 Mittelmann instances,
  of which only 7 pass the check. Its QP path is a barrier method
  whose solutions verify cleanly. We used cuOpt 26.8.0 with default settings apart
  from the tolerance and time limit.
* **PDLP (OR-Tools).** The same L2-relative termination rule applies, with a
  milder effect: 26 of its 306 "optimal" LP returns at :math:`10^{-4}` have
  infinity-norm residuals between :math:`10^{-3}` and :math:`10^{-2}`.
* **HiGHS.** Its QP solver is an active-set method and is slow on the larger
  QPs; its LP simplex is the fastest LP code in the comparison.
* **ProxQP.** Run with its default dense/sparse backend selection; it times
  out on many of the larger problems.
* **Clarabel.** Ran out of memory (64 GB) on ``equalG11``, which has a PSD
  block of order 801, so the other SDPLIB instances with PSD blocks of order
  500 or more were not attempted; all of them are counted as failures for it.
* **SDPA and CVXOPT.** Interior-point SDP solvers. SDPA does not take a time
  limit; its longest solve was 154 s, well inside the 900 s limit.
