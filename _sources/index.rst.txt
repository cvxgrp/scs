.. scs documentation master file, created by
   sphinx-quickstart on Sat Jul 24 12:54:37 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/scs_logo.png
  :width: 400
  :alt: SCS
  :align: center

.. title:: SCS

.. centered::
  **A fast, reliable, and open-source convex optimization solver.**

SCS (Splitting Conic Solver) is a numerical optimization package for solving
large-scale quadratic cone problems. It solves primal problems of the form

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

   * - :math:`s \in \mathbf{R}^m`
     - slack variable 

   * - :math:`y \in \mathbf{R}^m`
     - dual variable 

with data

.. list-table::
   :widths: 20 50
   :header-rows: 0

   * - :math:`A \in \mathbf{R}^{m \times n}`
     - sparse data matrix
   * - :math:`P \in \mathbf{S}_+^{n}`
     - sparse **symmetric positive semidefinite** matrix
   * - :math:`c \in \mathbf{R}^n`
     - dense primal cost vector 
   * - :math:`b \in \mathbf{R}^m`
     - dense dual cost vector 
   * - :math:`\mathcal{K} \subseteq \mathbf{R}^m`
     - nonempty, closed, convex cone, see :ref:`cones`
   * - :math:`\mathcal{K}^* \subseteq \mathbf{R}^m`
     - dual cone to :math:`\mathcal{K}`

At termination SCS will either return points :math:`(x^\star,y^\star,s^\star)` that satisfies
the :ref:`optimality conditions <solver>` to the desired accuracy, or a certificate
of :ref:`primal or dual infeasibility` to the designated infeasibility accuracy.

The current version is |version|

The code is freely available on `GitHub <https://github.com/cvxgrp/scs>`_. 

GPU

Features
________

.. glossary::

  Fast 
    SCS only requires a single sparse matrix factorization    

  Flexible 
    Cones

  Free and open source
    SCS is distributed under the permissive `MIT license <https://github.com/cvxgrp/scs/blob/master/LICENSE.txt>`_.

  Accelerated
    SCS includes Anderson acceleration that can accelerate convergence to high accuracy solutions.
 
  Detects infeasibility 
    SCS uses a homogeneous embedding and so can robustly and reliably detect infeasible problems.

  Language bindings
    SCS has bindings to many languages, including C, Python, Julia, R, MATLAB, and Ruby.
   
  Warm starts
    It can be easily warm-started and the matrix factorization can be cached.

  Matrix-free
    SCS has an indirect linear system solver

  Supported 
    SCS is a supported solver in parser-solvers CVX, CVXPY, YALMIP, 

  Battle-tested
    SCS was the first open-source solver available based on ADMM and is in wide-usage in industry and academia.



**Development.**

SCS is a community project, built from the contributions of many
researchers and engineers. The primary maintainer is
`Brendan O'Donoghue <https://bodono.github.io/>`_.
We appreciate all contributions. To get involved, see our :doc:`contributing
guide </contributing/index>`.


.. toctree::
   :hidden:
   :maxdepth: 2

   solver/index
   install/index
   help/index
   api/index
   linear_solver/index
   blas_lapack/index
   examples/index
   contributing/index
   citing/index

