:orphan:

.. This is the master document of the PDF user guide only (see latex_documents
   in conf.py). It is not linked from the site, whose navigation is the root
   toctree in index.rst. Pages are listed here selectively: the guide covers
   the C and Python interfaces and the material every user needs, and points
   to the site for the other language interfaces and the worked examples.

SCS User Guide
==============

This guide describes how to use SCS, the Splitting Conic Solver, through its
C and Python interfaces: how to install it, how to pose a problem in the form
it accepts, what its settings do, how to choose a linear-solver backend, and
what to do when a solve does not go as expected. The algorithm itself is
described near the end, as background. The release number on the title page
identifies the version of SCS the guide describes.

The MATLAB, Julia, R, Ruby and JavaScript interfaces, and worked examples
with their solver output, are documented on the website,
https://www.cvxgrp.org/scs/, from which this guide is generated.

To cite SCS, please cite the published papers listed in the final chapter
rather than this guide.

.. toctree::
   :maxdepth: 2

   install
   api
   /linear_solver/index
   /blas_lapack/index
   best_practices
   /help/index
   /algorithm/index
   /citing/index
