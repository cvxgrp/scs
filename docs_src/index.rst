.. scs documentation master file, created by
   sphinx-quickstart on Sat Jul 24 12:54:37 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/scs.png
  :width: 400
  :alt: SCS
  :align: center

Welcome to SCS
====================

**Convex optimization, for everyone.**
SCS (`splitting conic solver`) is a numerical optimization package for solving
large-scale convex cone problems, based on our paper.

.. math::
  \begin{array}{lcr}
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T P x + c^T x \\
    \mbox{subject to} &  A x + s = b\\
                      &  s \in \mathcal{K}
  \end{array}
  &\quad&
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T P x + c^T x \\
    \mbox{subject to} &  A x + s = b\\
                      &  s \in \mathcal{K}
  \end{array}
  \end{array}

**Current version**
The current version is |version|

**Code**
TODO

**Features**
TODO

**License**
TODO


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
   examples/index
   contributing/index
   citing/index

