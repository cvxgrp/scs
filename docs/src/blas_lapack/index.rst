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

Parallel projection of PSD blocks
---------------------------------
By default SCS projects the PSD blocks of a problem one after another, and a
multithreaded BLAS is free to use its threads inside each block. That is the
case for every build, with or without OpenMP, that is not compiled with
:code:`USE_PARALLEL_PSD`.

Compiled with :code:`USE_OPENMP=1 USE_PARALLEL_PSD=1`, SCS projects the real
and complex PSD blocks of a problem in parallel, each thread calling LAPACK for
its own blocks. For a problem with many blocks this can be several times
faster, but only if the BLAS then runs on one thread inside each thread of SCS.
That depends on the BLAS, which is why this is an option for whoever builds SCS
and chooses the BLAS, and not the default:

* A BLAS on the same OpenMP runtime as SCS runs on one thread there by itself
  (OpenMP builds of OpenBLAS and BLIS, MKL on the matching threading layer),
  because OpenMP does not nest parallel regions by default.
* A BLAS with threads of its own (a pthreads build of OpenBLAS), or on another
  OpenMP runtime than SCS, keeps its threads. The threads of the two then
  multiply, and the projection can be slower than without the option, unless
  the BLAS is set to one thread for the whole process
  (:code:`OPENBLAS_NUM_THREADS=1`, :code:`MKL_NUM_THREADS=1`). A pthreads
  OpenBLAS on one thread still takes a lock in its calls, which limits the
  gain.
* Accelerate keeps its own threads as well. On Apple silicon the parallel
  projection was nevertheless measured to be faster than the serial one with
  Accelerate at its default settings.

SCS guards against two of these cases where it can (see below), but it cannot
see every BLAS.

Each block goes through the same LAPACK calls as in the serial projection, so
with single-threaded BLAS the result is the same bit for bit, whatever the
number of threads (unless the rounding of the BLAS depends on the memory
alignment of its arguments, as it can with MKL).

The number of threads is fixed when the solver is set up. It is the OpenMP
thread count at that point (as set by :code:`OMP_NUM_THREADS` or
:code:`omp_set_num_threads`), at most the number of blocks larger than
:math:`1 \times 1`. A solve uses fewer threads if the OpenMP thread count has
been lowered since the setup. With :code:`verbose` set, the header gives the
number of threads of a parallel projection.

The blocks are still projected one after another when

* the total work is small, or one block holds more than half of it (the
  parallel projection would then gain less than a factor of two, and a
  multithreaded BLAS may do better on that block),
* the BLAS is a pthreads build of OpenBLAS set to more than one thread, where
  SCS can see that (see below), or
* the solver is set up inside an active OpenMP parallel region, where the
  threads of the caller are the parallelism.

What SCS does about the threads of the BLAS:

* MKL on another OpenMP runtime (for example its default Intel threading layer,
  when SCS runs on the GNU runtime): SCS sets the thread-local thread count of
  MKL to one during the projection. This goes through a weak reference to
  :code:`MKL_Set_Num_Threads_Local`, which needs an ELF platform such as Linux
  and the MKL symbols visible to SCS. Where that is not the case (Windows,
  macOS, a static MKL, MKL loaded privately by another library), set
  :code:`MKL_NUM_THREADS=1`, or select the threading layer that matches the
  runtime of SCS (:code:`MKL_THREADING_LAYER=GNU`).
* A pthreads build of OpenBLAS with more than one thread: on ELF platforms SCS
  checks for it at setup and then projects the blocks one after another. Set
  :code:`OPENBLAS_NUM_THREADS=1` to get the parallel projection. The check
  needs the functions :code:`openblas_get_parallel` and
  :code:`openblas_get_num_threads` visible to SCS: it is off for a static
  :code:`libopenblas.a` (the weak references can stay unresolved), for an
  OpenBLAS with suffixed symbol names, and for one loaded privately by another
  library (libblastrampoline in Julia).
* In all other cases (other platforms, other libraries, another OpenMP
  runtime) SCS does nothing: the BLAS has to be on one thread, or the build
  left without :code:`USE_PARALLEL_PSD`.

Several threads call BLAS and LAPACK at once, so the library must be safe to
call that way. A single-threaded OpenBLAS (built with :code:`USE_THREAD=0`) is
safe only if it was built with :code:`USE_LOCKING=1` (see OpenBLAS's
:code:`Makefile.rule`). SCS does not check this; OpenBLAS 0.3.28 and later list
:code:`USE_LOCKING` in :code:`openblas_get_config()`. OpenBLAS also has a table
of memory buffers of a size fixed when it is built (by :code:`NUM_THREADS`, and
at least 50), and a thread holds one during a BLAS call: with more threads than
buffers, OpenBLAS 0.3.18 and later print a warning and allocate more, and
earlier versions fail. The OpenMP thread count bounds the number of threads.

A parallel projection takes more memory: each thread after the first has its
own eigendecomposition workspace, with matrix buffers for the largest block
that thread projects and LAPACK work arrays of the serial sizes. The extra
memory is therefore at most the PSD workspace of the serial projection times
the number of threads less one. That workspace is about :math:`2n^2` floats for
the largest real block of size :math:`n` and :math:`2n^2` complex numbers for
the largest complex block, plus LAPACK work arrays linear in :math:`n`.
