.. _matrices:

Data matrices
=============

* Both matrices :math:`A` and :math:`P` must be in `Compressed Sparse Column (CSC) format <https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html>`_ using zero-based indexing.  

* The order of the rows of :math:`A` must be in the order that the cones appear in the :ref:`cones table <cones>`.  

* The rows of :math:`A` corresponding to semidefinite cones (as well as the rows of :math:`b`) must be handled carefully, as described in :ref:`sdcone`.

* The matrix :math:`P` must be symmetric positive semidefinite and only the **upper triangular** part of :math:`P` should be passed in.

