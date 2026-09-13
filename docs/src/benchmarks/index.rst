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
interior-point codes and clearly ahead of the other first-order solvers, and
the cuDSS GPU backend is the fastest solver we tested on the largest quarter of
the QP problems. On the SDP test sets the interior-point solvers are faster on
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

**Tolerances.** Every solver was run at a requested relative tolerance of
:math:`10^{-4}` (the headline plots) and :math:`10^{-6}` (the high-accuracy
plots), passing the tolerance through each solver's own absolute and relative
settings. For SCS this means ``eps_abs = eps_rel = tol`` with the iteration limit
raised so that only the time limit can stop it; all other settings were left at
their defaults.

**Timing.** The reported time is wall-clock time for the solver call, including
any presolve, factorization and GPU transfer, but excluding reading the problem
from disk. QP and LP solves were limited to 300 s and SDP solves to 900 s. A
performance profile shows, for each solver, the fraction of problems solved
within a factor :math:`\tau` of the fastest solver on that problem; failures
never count as solved. The shifted geometric mean uses a shift of 1 s and
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

**Problem sets.** QP: Maros-Meszaros (138), QPLIB continuous convex subset, and
the ``qpbenchmark`` MPC set. LP: Netlib (feasible), Kennington, the MIPLIB 2017
LP relaxations up to 20 MB, and the Mittelmann LP set. SDP: SDPLIB and the
Mittelmann SDP set. "Largest quartile" means the quarter of each family with
the most nonzeros in the constraint matrix (plus the Hessian for QPs).

**Exclusions.** Clarabel cannot form the dense scaling block for PSD cones of
order 500 or more; those instances are counted as failures for it. SDPA
ignores time limits and was allowed to run to completion. The SDPLIB archive's
``maxG55`` and ``maxG60`` files are corrupt and were dropped for all solvers.
cuOpt's QP path is an interior-point method whose factorization failed with a
numerical error on a subset of the Maros-Meszaros problems; those count as
failures.

.. _bench_qp:

Quadratic programs
------------------

221 problems: Maros-Meszaros (138), the convex continuous QPLIB instances and
the ``qpbenchmark`` MPC set. On the largest quarter of these, SCS with cuDSS
has the lowest shifted geometric mean solve time of any solver, and SCS on the
CPU is within 10% of PIQP, the fastest interior-point code, with the same
number of verified solves. Over all 221 problems the interior-point solvers
PIQP and Clarabel are faster on the small instances, where an SCS solve is
dominated by fixed setup cost, but SCS solves as many problems as they do.

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

.. list-table:: QP: verified solves and shifted geometric mean time (s); all 221 problems / largest quartile (56)
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
   * - PIQP
     - 216
     - 2.1
     - 52
     - 9.7
     - 214
     - 2.5
     - 50
     - 11.4
   * - Clarabel
     - 213
     - 2.6
     - 49
     - 13.0
     - 212
     - 3.0
     - 49
     - 13.9
   * - SCS (GPU, cuDSS)
     - 210
     - 4.4
     - 53
     - 8.2
     - 204
     - 7.9
     - 47
     - 23.6
   * - SCS (CPU, MKL Pardiso)
     - 208
     - 5.0
     - 53
     - 10.5
     - 202
     - 7.3
     - 47
     - 25.0
   * - OSQP
     - 202
     - 7.0
     - 46
     - 23.4
     - 185
     - 15.9
     - 43
     - 45.6
   * - cuOpt (GPU)
     - 169
     - 22.1
     - 35
     - 58.8
     - 163
     - 25.2
     - 30
     - 80.8
   * - HiGHS
     - 162
     - 27.7
     - 24
     - 164.8
     - 144
     - 43.2
     - 17
     - 286.9
   * - ProxQP
     - 129
     - 69.8
     - 24
     - 212.2
     - 125
     - 76.6
     - 15
     - 425.8


.. _bench_lp:

Linear programs
---------------

345 problems: Netlib, Kennington, the MIPLIB 2017 LP relaxations up to 20 MB
and the Mittelmann LP set. SCS is not an LP solver and is not marketed as one,
so this comparison is included mainly to show that it is not a bad one: on the
largest quarter of the LP set SCS with cuDSS verifies the most solutions of any
solver and is second only to HiGHS (dual simplex) in geometric mean time,
ahead of the interior-point solvers PIQP and Clarabel. Both PDLP
implementations, which are first-order LP methods, trail SCS by a wide margin
under independent verification (see the notes below).

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

.. list-table:: LP: verified solves and shifted geometric mean time (s); all 345 problems / largest quartile (87)
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
     - 329
     - 5.3
     - 74
     - 27.1
     - 331
     - 4.8
     - 75
     - 25.0
   * - PIQP
     - 318
     - 8.5
     - 68
     - 40.6
     - 310
     - 9.8
     - 65
     - 43.0
   * - SCS (GPU, cuDSS)
     - 316
     - 10.8
     - 75
     - 33.9
     - 299
     - 16.5
     - 67
     - 49.0
   * - SCS (CPU, MKL Pardiso)
     - 315
     - 11.0
     - 72
     - 48.3
     - 295
     - 18.9
     - 62
     - 79.2
   * - Clarabel
     - 315
     - 11.2
     - 67
     - 65.3
     - 312
     - 12.4
     - 66
     - 72.7
   * - PDLP (OR-Tools)
     - 281
     - 23.9
     - 48
     - 140.5
     - 173
     - 129.9
     - 40
     - 227.8
   * - OSQP
     - 257
     - 41.0
     - 45
     - 205.9
     - 177
     - 132.0
     - 23
     - 484.3
   * - cuOpt (GPU)
     - 133
     - 164.5
     - 25
     - 269.2
     - 123
     - 188.3
     - 20
     - 346.3


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
     - 80
     - 21.7
     - 27
     - 12.3
     - 66
     - 53.5
     - 23
     - 35.5
   * - CVXOPT
     - 79
     - 38.4
     - 25
     - 35.3
     - 59
     - 108.0
     - 20
     - 88.4
   * - Clarabel
     - 48
     - 79.7
     - 13
     - 118.1
     - 47
     - 81.3
     - 11
     - 147.1
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


Notes on individual solvers
---------------------------

* **SCS.** A few SCS "solved" returns fail the independent residual check
  and count as failures above. Two of them are real false certificates
  (``QGROW7`` and ``QGROW22`` from Maros-Meszaros, and the MIPLIB relaxation
  ``neos-4413714-turia``): on these badly scaled problems the diagonal
  rescaling drives the internal scale to its floor, and the scaled residuals
  SCS monitors no longer track the unscaled ones. The rest (``PRIMALC8``,
  ``QSCFXM2``, ``greenbea``, the ``pilot`` family) are marginal, with
  residuals between one and four times the :math:`10^{-3}` cut. The same
  check promotes SCS solves that hit the time limit with residuals inside ten
  times the tolerance; those are counted as solved at the time limit, which
  is what makes the SCS curves on the SDP profile reach the right edge.
* **cuOpt.** cuOpt's LP path is PDLP with an L2-relative stopping rule, and on
  problems with many variable bounds (Kennington, MIPLIB) its returned duals
  often have large infinity-norm stationarity residuals even though cuOpt
  reports the solve as optimal, and even though its own reported absolute
  dual residual is in the tens or hundreds. Under the uniform check used here
  those solves count as failures; by its own status cuOpt reports 233 of the
  345 LPs optimal at :math:`10^{-4}`. Its QP path is a barrier method whose
  solutions verify cleanly. We used cuOpt 26.8.0 with default settings apart
  from the tolerance and time limit.
* **PDLP (OR-Tools).** The same L2-relative termination rule applies, with a
  milder effect: 26 of its 307 "optimal" LP returns at :math:`10^{-4}` have
  infinity-norm residuals between :math:`10^{-3}` and :math:`10^{-2}`.
* **HiGHS.** Its QP solver is an active-set method and is slow on the larger
  QPs; its LP simplex is the fastest LP code in the comparison.
* **ProxQP.** Run with its default dense/sparse backend selection; it times
  out on many of the larger problems.
* **Clarabel.** Cannot attempt PSD cones of order 500 or more (it forms a
  dense scaling block); those instances are counted as failures.
* **SDPA and CVXOPT.** Interior-point SDP solvers; SDPA has no time limit and
  ran to completion on every instance.
