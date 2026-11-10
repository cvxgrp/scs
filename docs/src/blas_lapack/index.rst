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
or complex PSD blocks, the blocks are projected in parallel. Each block goes
through the same LAPACK calls as in the serial projection, so with
single-threaded BLAS the result is the same bit for bit, whatever the number of
threads (unless the rounding of the BLAS depends on the memory alignment of its
arguments, as it can with MKL).

The number of threads is fixed when the solver is set up. It is the OpenMP
thread count at that point (as set by :code:`OMP_NUM_THREADS` or
:code:`omp_set_num_threads`), at most the number of blocks larger than
:math:`1 \times 1`, and at most :code:`SCS_PSD_NUM_THREADS` if that environment
variable is set to a positive number: :code:`SCS_PSD_NUM_THREADS=1` keeps the
projection serial. A solve uses fewer threads if the OpenMP thread count has
been lowered since the setup. With :code:`verbose` set, the header gives the
number of threads of a parallel projection.

The blocks are projected serially, as they are without OpenMP, when

* the total work is small, or one block holds more than half of it (the
  parallel projection would then gain less than a factor of two, and a
  multithreaded BLAS may do better on that block),
* the BLAS is a pthreads build of OpenBLAS set to more than one thread (see
  below), or
* the solver is set up inside an active OpenMP parallel region, where the
  threads of the caller are the parallelism.

In the parallel projection each thread calls LAPACK itself, so BLAS should not
start threads of its own there:

* OpenMP builds of OpenBLAS and BLIS, and MKL, run single threaded inside the
  parallel region when they use the same OpenMP runtime as SCS.
* MKL on another OpenMP runtime (for example its default Intel threading layer,
  when SCS runs on the GNU runtime) would start its own threads in each thread
  of SCS. SCS therefore sets the thread-local thread count of MKL to one during
  the projection. This goes through a weak reference to
  :code:`MKL_Set_Num_Threads_Local`, which needs an ELF platform such as Linux
  and the MKL symbols visible to SCS. Where that is not the case (Windows, a
  static MKL, MKL loaded privately by another library), set
  :code:`MKL_NUM_THREADS=1`, or select the threading layer that matches the
  runtime of SCS (:code:`MKL_THREADING_LAYER=GNU`).
* A pthreads build of OpenBLAS has no thread-local setting, and with more than
  one thread it would oversubscribe the cores. On ELF platforms SCS checks for
  one at setup and then projects serially: set :code:`OPENBLAS_NUM_THREADS=1`
  to get the parallel projection. The check needs the functions
  :code:`openblas_get_parallel` and :code:`openblas_get_num_threads` visible to
  SCS: it is off for a static :code:`libopenblas.a` (the weak references can
  stay unresolved), for an OpenBLAS with suffixed symbol names, and for one
  loaded privately by another library (libblastrampoline in Julia).
* In all other cases (other platforms, other libraries such as Accelerate), SCS
  assumes that BLAS is single threaded.

Where SCS cannot see or limit the threads of a multithreaded BLAS, set the BLAS
to one thread, or set :code:`SCS_PSD_NUM_THREADS=1`.

Several threads then call BLAS and LAPACK at once, so the library must be safe
to call that way. A single-threaded OpenBLAS (built with :code:`USE_THREAD=0`)
is safe only if it was built with :code:`USE_LOCKING=1` (see OpenBLAS's
:code:`Makefile.rule`). SCS does not check this; OpenBLAS 0.3.28 and later list
:code:`USE_LOCKING` in :code:`openblas_get_config()`. OpenBLAS also has a table
of memory buffers of a size fixed when it is built (by :code:`NUM_THREADS`, and
at least 50), and a thread holds one during a BLAS call: with more threads than
buffers, OpenBLAS 0.3.18 and later print a warning and allocate more, and
earlier versions fail. :code:`SCS_PSD_NUM_THREADS` bounds the number of threads.

A parallel projection takes more memory: each thread after the first has its
own eigendecomposition workspace, with matrix buffers for the largest block
that thread projects and LAPACK work arrays of the serial sizes. The extra
memory is therefore at most the PSD workspace of the serial projection times
the number of threads less one. That workspace is about :math:`2n^2` floats for
the largest real block of size :math:`n` and :math:`2n^2` complex numbers for
the largest complex block, plus LAPACK work arrays linear in :math:`n`.
