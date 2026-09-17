.. scs documentation master file, created by
   sphinx-quickstart on Sat Jul 24 12:54:37 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/scs_logo_transparent.png
  :width: 400
  :alt: SCS
  :align: center
  :target: https://github.com/cvxgrp/scs

.. title:: SCS

.. centered::
  **A fast, reliable, and open-source convex cone solver.**

SCS (Splitting Conic Solver) is a numerical optimization package for solving
large-scale convex quadratic cone problems. The code is freely available on
`GitHub <https://github.com/cvxgrp/scs>`_.  It solves primal-dual problems of
the form

.. math::
  \begin{array}{lcr}
  \begin{array}{ll}
  \mbox{minimize} & (1/2)x^\top  P x + c^\top  x\\
  \mbox{subject to} &  Ax + s = b\\
    & s \in \mathcal{K}
  \end{array}
  &\quad&
  \begin{array}{ll}
  \mbox{maximize} & -(1/2)x^\top  P x - b^\top  y\\
  \mbox{subject to} &  Px + A^\top y + c = 0\\
    & y \in \mathcal{K}^*
  \end{array}
  \end{array}

over variables

.. list-table::
   :widths: 20 20
   :header-rows: 0

   * - :math:`x \in \mathbf{R}^n`
     - primal variable

   * - :math:`y \in \mathbf{R}^m`
     - dual variable

   * - :math:`s \in \mathbf{R}^m`
     - slack variable

with data

.. list-table::
   :widths: 20 50
   :header-rows: 0

   * - :math:`A \in \mathbf{R}^{m \times n}`
     - sparse data matrix, see :ref:`matrices`
   * - :math:`P \in \mathbf{S}_+^{n}`
     - sparse, symmetric positive semidefinite matrix
   * - :math:`c \in \mathbf{R}^n`
     - dense primal cost vector
   * - :math:`b \in \mathbf{R}^m`
     - dense dual cost vector
   * - :math:`\mathcal{K} \subseteq \mathbf{R}^m`
     - nonempty, closed, convex cone, see :ref:`cones`
   * - :math:`\mathcal{K}^* \subseteq \mathbf{R}^m`
     - dual cone to :math:`\mathcal{K}`

At termination SCS will either return points :math:`(x^\star,y^\star,s^\star)` that satisfies
the :ref:`optimality conditions <optimality>` to the desired accuracy, or a certificate
of :ref:`primal or dual infeasibility <infeasibility>` to the designated infeasibility accuracy.

Features
--------

* **Efficient**: Designed to scale to large problems.
* **Flexible**: Supports quadratic objectives and a large range of :ref:`cones <cones>`.
* **Free and open source**: Distributed under the permissive `MIT license <https://github.com/cvxgrp/scs/blob/master/LICENSE.txt>`_.
* **Detects infeasibility**: Robustly and reliably detects :ref:`infeasible <infeasibility>` problems.
* **Interfaces**: Bindings for many :ref:`languages <interfaces>`, including C, Python, Julia, R, MATLAB, Ruby, and JavaScript via WebAssembly.
* **Warm starts**: Easily :ref:`warm-started <warm_start>`, and the matrix factorization can be cached.
* **Matrix-free**: Optionally use an :ref:`indirect linear system solver <indirect>`, or a :ref:`GPU version <gpu_indirect>`.
* **Supported**: A supported solver in `CVX <http://cvxr.com/cvx/>`_, `CVXPY <https://github.com/cvxgrp/cvxpy>`_, `YALMIP <https://github.com/johanlofberg/YALMIP>`_, `Convex.jl <https://github.com/jump-dev/Convex.jl>`_  and `JuMP <https://github.com/jump-dev/JuMP.jl>`_.
* **Accelerated**: Includes :ref:`acceleration <acceleration>` that can improve convergence to high accuracy.
* **Battle-tested**: The first ADMM-based solver available, and in wide usage.

Performance
-----------

SCS 3.3 is built for large problems. On the Mittelmann set of large LPs it
solves more instances than any other open-source solver; with the cuDSS GPU
backend it is the fastest solver we have measured on the largest QPs and
second only to HiGHS on the largest LPs; and on the CPU it matches the best
interior-point solvers problem for problem at a relative tolerance of
:math:`10^{-4}`, while also handling the second-order, semidefinite,
exponential and power cones that QP solvers cannot. Every solution is verified
independently against the same residual test, so solvers with different ideas
of "tolerance" are compared fairly. SCS also detects when there is no
solution, returning certificates of infeasibility or unboundedness that are
verified the same way. The full methodology, all test sets including SDPs
and the infeasible problems, and the raw data are on the :ref:`benchmarks
page <benchmarks>`.

.. figure:: files/bench/landing_grid.png
   :width: 100 %
   :align: center
   :alt: Performance profiles and geometric mean solve times on the largest QPs, largest LPs and the Mittelmann LP set

   Largest quarter of the Maros-Meszaros and QPLIB QPs (top), largest quarter
   of the Kennington and MIPLIB-relaxation LPs (middle) and the Mittelmann
   large-LP set (bottom), at tolerance
   :math:`10^{-4}`. Left: fraction of problems solved within a factor
   :math:`\tau` of the fastest solver. Right: shifted geometric mean solve
   time with failures charged three times the time limit, and the number of
   verified solves.

.. figure:: files/bench/infeas_pair.png
   :width: 100 %
   :align: center
   :alt: Performance profile and geometric mean time to a verified certificate on the 29 infeasible Netlib LPs

   Infeasibility detection on the 29 infeasible Netlib LPs, for the solvers
   that return certificates: a solve is a certificate of infeasibility or
   unboundedness verified from the problem data. These are small problems
   (a median of 460 variables), so this tests detection rather than speed
   at scale; settings are described on the :ref:`benchmarks page
   <bench_infeasible>`.


Development
-----------

SCS is a community project, built from the contributions of many
researchers and engineers. The primary maintainer is
`Brendan O'Donoghue <https://bodono.github.io/>`_.
We appreciate all contributions. To get involved, see our :doc:`contributing
guide </contributing/index>`.



.. toctree::
   :hidden:
   :maxdepth: 2

   algorithm/index
   benchmarks/index
   api/index
   install/index
   linear_solver/index
   blas_lapack/index
   examples/index
   contributing/index
   help/index
   citing/index

