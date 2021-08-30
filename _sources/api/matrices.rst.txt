.. _matrices:

Matrices
--------

The matrices :math:`A` and :math:`P` must be in `Compressed Sparse Column (CSC)
format <https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html>`_ using zero-based
indexing.  The order of the rows of :math:`A` must be in the order that the
cones appear in the table :ref:`here <cones>`.  The matrix :math:`P` must be
**symmetric positive semidefinite** and only the **upper triangular** part of
:math:`P` should be passed in.

