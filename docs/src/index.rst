.. scs documentation master file, created by
   sphinx-quickstart on Sat Jul 24 12:54:37 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SCS 3.0.0
====================

**Convex optimization, for everyone.**
SCS (`splitting conic solver`) is a numerical optimization package for solving
large-scale convex cone problems, based on our paper.

.. image:: _static/scs.png
  :width: 400
  :alt: SCS
  :align: center

.. math::
     \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T P x + q^T x \\
    \mbox{subject to} & l \leq A x \leq u
  \end{array}


**Development.**

SCS is a community project, built from the contributions of many
researchers and engineers.

SCS is developed and maintained by
`Brendan O'Donoghue <https://bodono.github.io/>`_ with many others
contributing significantly. 
We appreciate all contributions. To get involved, see our :doc:`contributing
guide </contributing/index>`.

.. toctree::
      :hidden:

   install/index

.. toctree::
       :maxdepth: 3
    :hidden:

    tutorial/index

.. toctree::
      :hidden:

   examples/index

.. toctree::
      :hidden:

   API Documentation <api_reference/cvxpy>

.. toctree::
      :maxdepth: 1
   :hidden:

   faq/index

.. toctree::
      :hidden:

   citing/index

.. toctree::
      :hidden:

   contributing/index

.. toctree::
      :hidden:

   related_projects/index

.. toctree::
      :hidden:

   updates/index

.. toctree::
      :hidden:

   short_course/index

.. toctree::
      :hidden:

   license/index
