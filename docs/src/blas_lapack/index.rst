.. _blas_lapack:

BLAS and LAPACK
===============
BLAS and LAPACK are dense numerical linear algebra packages. SCS uses these
libraries in two places:

* To compute the eigen-decomposition for the semidefinite :ref:`cone <cones>` projection
* To solve the linear system in :ref:`Anderson acceleration <acceleration>`

Therefore compiling with BLAS / LAPACK **is optional**. If you are
not interested in solving SDPs or using acceleration then there is no need to
use these libraries. To compile without these libraries you can set the
:ref:`compiler flag <compile_flags>` :code:`USE_LAPACK` to :code:`0`, e.g., if
:ref:`installing <install>` using :code:`make`:

.. code:: bash

  make USE_LAPACK=0 

If you do want to solve SDPs or use Anderson acceleration, then you will need
to install BLAS and LAPACK libraries (these are pre-installed in most machines).
If calling SCS from another language (Python, MATLAB etc.) then these libraries
should be pre-installed and SCS will try to link against them. Otherwise
you may need to install a copy yourself. A good library to start with is
`OpenBLAS <https://www.openblas.net/>`_, which contains both BLAS and LAPACK.

There are many different BLAS and LAPACK libraries that conform to the same API.
Finding one that is optimized for your machine can make a big difference to
the speed of the operations in practice. If the speed of the SDP projection
or the acceleration step is a bottleneck you can experiment with faster
libraries like `MKL <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html>`_ 
or `ATLAS <http://math-atlas.sourceforge.net/>`_.


