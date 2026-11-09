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



.. _blas_lapack_openmp:

OpenMP and BLAS threads
-----------------------
When SCS is compiled with :code:`USE_OPENMP=1` and a problem has several real
or complex PSD blocks, the blocks are projected in parallel over
:code:`OMP_NUM_THREADS` threads (read when the solver is set up); with
single-threaded BLAS the result is the same, bit for bit, as the serial
projection. Each thread calls LAPACK itself, so BLAS should not start threads of
its own: set :code:`OPENBLAS_NUM_THREADS=1` for a pthreads build of OpenBLAS
(on Linux and other ELF platforms SCS detects one and then projects serially;
elsewhere, and for other BLAS libraries, it assumes BLAS is single threaded).
The check needs OpenBLAS linked dynamically: with a static :code:`libopenblas.a`
the weak references can stay unresolved, and the check is off.

Several threads then call BLAS and LAPACK at once, so the library must be safe
to call that way. A single-threaded OpenBLAS (built with :code:`USE_THREAD=0`)
is safe only if it was built with :code:`USE_LOCKING=1` (see OpenBLAS's
:code:`Makefile.rule`); SCS cannot detect this at run time.

The blocks are projected serially, with no extra memory, when the total work is
small or one block holds more than half of it. Otherwise each extra thread has
its own eigendecomposition workspace: matrix buffers for the largest block that
thread projects, plus LAPACK work arrays of the serial sizes. The extra memory is
therefore at most :code:`OMP_NUM_THREADS - 1` times the PSD workspace of the
serial projection, which is about :math:`2n^2` floats for the largest real block
of size :math:`n` and :math:`2n^2` complex numbers for the largest complex block,
plus LAPACK work arrays linear in :math:`n`.
