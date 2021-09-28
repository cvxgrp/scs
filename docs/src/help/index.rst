.. _help:

Help
====

Currently the easiest way to get support is to file a `GitHub issue
<https://github.com/cvxgrp/scs/issues>`_ or to `email us
<mailto:splitting.conic.solver@gmail.com>`_ (which may be slower).

If you are have a problem that SCS struggles to solve you can set the
:code:`write_data_filename` field in the :ref:`settings <settings>` and SCS will
dump a file containing the problem data to disk under that filename. Zip the
file and `email it to us <mailto:splitting.conic.solver@gmail.com>`_ or attach
it to a GitHub issue. This makes it much easier for us to reproduce the problem.

A common cause of issues is not linking :ref:`BLAS/LAPACK libraries
<blas_lapack>` correctly. If you are having this issue please search for
resources on installing and linking these libraries first. You can try `OpenBLAS
<https://www.openblas.net/>`_ if you need a BLAS library.

