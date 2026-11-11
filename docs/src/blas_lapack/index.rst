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
A problem with several positive semidefinite (PSD) blocks spends much of its
time projecting them, and SCS can project several blocks at the same time on
different threads. It does so **only in a build that asks for it**, with the
build option :code:`USE_PARALLEL_PSD`. Without the option the blocks are
projected one after another. That is the default, and it is how the Python
wheels and the Julia packages are built.

The option is off by default because it is faster only with a BLAS that is set
up for it, and slower with one that is not. This section explains when it
helps, how to turn it on, and how to set up the BLAS.

When it helps
~~~~~~~~~~~~~
* It helps with **several PSD blocks of comparable size**, many small ones or a
  few large ones. On 4 threads, 4 to 16 blocks of size 256 to 1024 were
  projected about 3 times faster, and a thousand blocks of size 10 to 32 about
  3.5 times faster.
* It does not help with **one block, or one block much larger than the rest**:
  a parallel projection is no faster than its largest block. SCS recognizes
  such a problem and projects its blocks one after another, as it does without
  the option.
* It speeds up the projection, not the rest of a solve. If most of the time
  goes into the linear system, the solve as a whole gains little.

How to turn it on
~~~~~~~~~~~~~~~~~
Build SCS with OpenMP and with the option:

.. code:: bash

  make USE_OPENMP=1 USE_PARALLEL_PSD=1

or, with CMake,

.. code:: bash

  cmake -DUSE_OPENMP=ON -DUSE_PARALLEL_PSD=ON ..

Then set the OpenMP thread count to the number of cores to use:

.. code:: bash

  export OMP_NUM_THREADS=4

On a machine with fast and slow cores, count only the fast ones. Every thread
gets an equal share of the blocks, so a thread on a slow core holds the others
up.

Setting up the BLAS
~~~~~~~~~~~~~~~~~~~
Threads can work on the projection at two levels. SCS can project several
blocks at once, and the BLAS can use several threads inside the
eigendecomposition of a single block. **Only one of the two should be active at
a time.** With the option on, every thread of SCS calls the BLAS, so the BLAS
has to run on one thread there. If it does not, each thread of SCS starts the
threads of the BLAS as well. There are then more threads than cores, and the
projection is slower than it should be: in one measurement the process had 16
threads where 4 were asked for, and ran at half the speed.

Whether the BLAS stays on one thread depends on where it gets its threads from:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - Your BLAS
     - Inside the threads of SCS
     - What to do
   * - Built with OpenMP, on the same OpenMP runtime as SCS: an OpenMP build of
       OpenBLAS or BLIS, NVPL, or MKL on the threading layer that matches the
       compiler of SCS
     - Runs on one thread by itself, because OpenMP does not nest parallel
       regions by default
     - Nothing. **This is the setup to use**
   * - MKL on Intel's OpenMP runtime (its default), with SCS compiled by GCC
     - Would start its own threads. On Linux, with MKL as a shared library,
       SCS holds it to one thread during the projection
     - Nothing on Linux. Elsewhere set :code:`MKL_THREADING_LAYER=GNU` or
       :code:`MKL_NUM_THREADS=1`
   * - A pthreads build of OpenBLAS, which has a thread pool of its own. The
       default :code:`libopenblas-dev` of Debian and Ubuntu is one
     - Keeps its threads. On Linux SCS sees that and projects the blocks one
       after another instead
     - Use an OpenMP build of OpenBLAS if you can. Otherwise set
       :code:`OPENBLAS_NUM_THREADS=1`: the projection is then parallel, but
       this OpenBLAS takes a lock in its calls and gains less
   * - Built with OpenMP, on a different OpenMP runtime from SCS (the BLAS on
       LLVM's runtime and SCS compiled by GCC, say)
     - Cannot see the threads of SCS and starts its own in each of them
     - Build SCS and the BLAS with the same compiler family. No variable puts
       only the BLAS on one thread: :code:`OMP_NUM_THREADS` is read by both
   * - Accelerate (macOS)
     - Keeps a pool of threads of its own, which stays almost idle
     - Nothing. At its default settings the parallel projection was 2.3 to 3.3
       times faster than the serial one on 4 threads

Putting a BLAS on one thread with an environment variable holds for the whole
process: that BLAS then runs on one thread wherever else it is used, in SCS
(the dense linear solver) and in your own code. A BLAS on the OpenMP runtime of
SCS has no such cost. The OpenMP thread count is then one budget, which SCS
uses during a parallel projection and the BLAS uses everywhere else.

Checking that it works
~~~~~~~~~~~~~~~~~~~~~~
* With :code:`verbose` set, the header of a solve says how many threads project
  the PSD blocks, for example ``psd cones projected in parallel by 4 threads``.
  Without that line the blocks are projected one after another; the next part
  lists the reasons.
* Count the threads of the process during a solve (``top -H``, or
  ``ps -o nlwp= -p <pid>`` on Linux). There should be about as many as the
  OpenMP thread count. Several times more means that the BLAS is running its
  own threads inside those of SCS.
* Compare the time of a solve with that of a build without the option. If it
  is not faster, leave the option off.

How SCS uses the threads
~~~~~~~~~~~~~~~~~~~~~~~~
The number of threads is fixed when the solver is set up. It is the OpenMP
thread count at that point (as set by :code:`OMP_NUM_THREADS` or
:code:`omp_set_num_threads`), at most the number of blocks larger than
:math:`1 \times 1`. It is also bounded by the amount of work, because each
thread takes a few microseconds to start: a projection of some tens of
microseconds gets two threads, and one of a millisecond or more gets all of
them. A solve uses fewer threads if the OpenMP thread count has been lowered
since the setup.

Even with the option, the blocks are projected one after another, and the BLAS
is free to use its own threads, when

* the threads would gain little: the blocks are too little work, or one block
  holds most of it (the parallel projection is no faster than its largest
  block),
* a block larger than :math:`512 \times 512` holds more than half of the work
  (the parallel projection would then gain less than a factor of two, and a
  multithreaded BLAS may do better on that block),
* the BLAS is a pthreads build of OpenBLAS set to more than one thread, where
  SCS can see that (see below), or
* the solver is set up inside an active OpenMP parallel region, where the
  threads of the caller are the parallelism.

What SCS checks, and where it cannot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCS guards against two of the cases in the table, on some platforms only. That
is the reason the option is left to whoever builds SCS and knows its BLAS.

* MKL on another OpenMP runtime: SCS sets the thread-local thread count of MKL
  to one during the projection. This goes through a weak reference to
  :code:`MKL_Set_Num_Threads_Local`, which needs an ELF platform such as Linux
  and the MKL symbols visible to SCS. It does not work on Windows or macOS,
  with a static MKL, or with an MKL loaded privately by another library.
* A pthreads build of OpenBLAS with more than one thread: on ELF platforms SCS
  checks for it at setup and then projects the blocks one after another. The
  check needs the functions :code:`openblas_get_parallel` and
  :code:`openblas_get_num_threads` visible to SCS: it is off for a static
  :code:`libopenblas.a` (the weak references can stay unresolved), for an
  OpenBLAS with suffixed symbol names, and for one loaded privately by another
  library (libblastrampoline in Julia).
* In all other cases (other platforms, other libraries, another OpenMP
  runtime) SCS does nothing and assumes that the BLAS runs on one thread.

Thread safety, memory and reproducibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

Each block goes through the same LAPACK calls as in the serial projection, so
with single-threaded BLAS the result is the same bit for bit, whatever the
number of threads (unless the rounding of the BLAS depends on the memory
alignment of its arguments, as it can with MKL).
