.. _blas_lapack:

BLAS and LAPACK
===============
BLAS and LAPACK are dense numerical linear algebra packages. SCS uses these
libraries in the following places:

* To compute the eigen-decomposition for the semidefinite :ref:`cone <cones>`
  projection, and the eigen-decompositions and singular value decompositions
  for the :ref:`spectral cones <spectral_cones>`
* To solve the linear system in :ref:`Anderson acceleration <acceleration>`
* To form and factorize the matrix in the :ref:`dense direct <dense>` linear
  solver
* To compute the eigenvectors that the :ref:`indirect <indirect>` linear
  solver uses for deflation
* For the basic vector operations (norms, inner products, scaling and adding
  vectors), which are plain C loops otherwise

Compiling with BLAS / LAPACK **is optional**: SCS builds and runs without
them. A build without them cannot solve SDPs (the setup fails for a
semidefinite cone larger than :math:`1 \times 1`), cannot include the spectral
cones or the :ref:`parallel PSD projection <blas_lapack_openmp>`, does not
apply Anderson acceleration (whatever the value of
:code:`acceleration_lookback`), does not have the dense direct linear solver,
and runs the indirect linear solver without deflation. If you need none of
these then there is no need to use these libraries. To compile without these
libraries you can set the
:ref:`compiler flag <compile_flags>` :code:`USE_LAPACK` to :code:`0`, e.g., if
:ref:`installing <install>` using :code:`make`:

.. code:: bash

  make USE_LAPACK=0 

If you do want any of the above, then you will need
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


.. _blas_lapack_threads:

Number of BLAS threads
----------------------
SCS projects each PSD block with a LAPACK eigendecomposition (:code:`syevr`,
or :code:`heevr` for a complex block) and a BLAS product (:code:`syrk` or
:code:`herk`). A multithreaded BLAS splits these calls over its own threads,
by default about one per core. Waking and synchronizing the threads has a
fixed cost, so blocks of size up to a few hundred gain nothing from them,
and with some libraries they are slower. Large blocks gain. The number of
threads is set by an environment variable of the BLAS, so you can choose.

What we measured
~~~~~~~~~~~~~~~~
Time of the PSD projection per iteration on the SDPLIB problems, relative to
the BLAS on one thread (geometric mean, range in brackets; below 1 is
faster). The machine is an NVIDIA GB10 with 10 fast and 10 slow ARM cores,
the libraries NVPL and a pthreads build of OpenBLAS 0.3.31. The runs on one
and on 4 threads were pinned to fast cores; 20 threads is the default of both
libraries there.

.. list-table::
   :header-rows: 1
   :widths: 24 19 19 19 19

   * - Size of the largest PSD block
     - NVPL, 4 threads
     - NVPL, 20 threads
     - OpenBLAS, 4 threads
     - OpenBLAS, 20 threads
   * - 37 to 130 (25 problems)
     - 1.34 [1.03, 1.93]
     - 2.27 [1.20, 3.83]
     - 1.02 [0.99, 1.05]
     - 1.18 [1.01, 1.90]
   * - 150 to 300 (17)
     - 1.07 [1.01, 1.15]
     - 2.45 [1.75, 6.34]
     - 0.94 [0.85, 1.01]
     - 1.55 [1.02, 2.10]
   * - 500 (8)
     - 0.90 [0.89, 0.93]
     - 1.32 [1.28, 1.37]
     - 0.76 [0.74, 0.79]
     - 1.05 [0.84, 1.23]
   * - 800 to 2000 (8)
     - 0.68 [0.50, 0.88]
     - 0.82 [0.49, 1.10]
     - 0.61 [0.45, 0.88]
     - 0.59 [0.33, 0.89]

* control7 (blocks of size 70 and 35) took 0.51 ms per projection with NVPL
  on one thread, 0.98 ms on 4 and 1.74 ms on 20, and a whole iteration
  1.29 ms on one thread against 2.63 ms on 20.
* The pthreads OpenBLAS lost little time on small blocks but kept cores busy
  for nothing: at its default, the 10 problems with blocks of size 37 to 124
  that we checked used 3 to 19 cores (310% to 1920% CPU), where one thread
  was as fast or faster. Those cores are taken from anything else on the
  machine, such as other solves.
* The 20-thread columns include the slow cores and depend on where the system
  placed the threads: problems with tiny blocks, which 4 threads did not
  change, were also up to 2 times slower there. So they show what the
  default can cost on such a machine, not the cost of the threads alone.

When to set one thread
~~~~~~~~~~~~~~~~~~~~~~
Put the BLAS on one thread for an SDP whose PSD blocks are all of size below a
few hundred (up to 200, 4 threads gained nothing with either library), in
particular when several solves run at the same time, since each would start
threads for every core. Keep the threads when

* the problem has PSD blocks of size 500 or more, which were faster on 4 threads
  with both libraries;
* the linear system is large and solved with MKL Pardiso, which runs on MKL's
  threads (see the threads paragraph of :ref:`MKL Pardiso <mkl>`).

A problem with several PSD blocks can instead use threads across the blocks,
with the BLAS on one thread inside each: see :ref:`blas_lapack_openmp`.

Which variable to set
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Your BLAS
     - Set
     - Notes
   * - OpenBLAS built with pthreads, its default (the default
       :code:`libopenblas-dev` of Debian and Ubuntu is one)
     - :code:`OPENBLAS_NUM_THREADS=1`
     - When that is unset it reads :code:`OMP_NUM_THREADS`, so
       :code:`OMP_NUM_THREADS=1` works too
   * - A BLAS built with OpenMP: an OpenMP build of OpenBLAS, NVPL
     - :code:`OMP_NUM_THREADS=1`
     - An OpenMP build of OpenBLAS ignores :code:`OPENBLAS_NUM_THREADS`
   * - MKL
     - :code:`MKL_NUM_THREADS=1`
     - Read before :code:`OMP_NUM_THREADS`. It also puts MKL Pardiso on one
       thread. A sequential MKL, such as the MKL backend of the x86-64 Linux
       Python wheels, always runs on one thread. Not measured in the
       projection
   * - Accelerate (macOS)
     - Nothing
     - :code:`VECLIB_MAXIMUM_THREADS=1` made no difference to the projection
       on an Apple M4 (82 SDPLIB problems: time ratios 0.89 to 1.11, 1.01 on
       average)

Set the variable before the program starts (in Python, at the latest before
importing NumPy or SCS): the libraries read it when they start. It holds for
every call to that library in the process, including SCS's dense direct
linear solver and your own code.

:code:`OMP_NUM_THREADS` is also SCS's own thread count in a build with OpenMP
(:code:`USE_OPENMP`; the pre-built Python wheels are built without it).
Setting it to 1 there also stops SCS's own threads, the :ref:`parallel PSD
projection <blas_lapack_openmp>` included, so use the library's own variable
where it has one.

To check, watch the CPU use of the process during a solve (``top``). With the
BLAS on one thread and SCS built without OpenMP it stays near 100%; several
hundred percent on a problem with small PSD blocks means that the BLAS threads
are running.


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

With CMake the option takes effect only in the libraries that CMake compiles
with OpenMP: :code:`scsdir`, :code:`scsindir`, :code:`scsdense`,
:code:`scsaccel` and :code:`scsgpuindir`. CMake compiles the MKL Pardiso and
cuDSS libraries (:code:`scsmkl` and :code:`scscudss`) without OpenMP, so they
project the blocks one after another whatever the option is set to. The
Makefile compiles every library with OpenMP, these two included.

With :code:`-DUSE_PARALLEL_PSD=ON`, if CMake does not find OpenMP the
configuration stops with an error, because no library would then have the
parallel projection. OpenMP that is installed where CMake does not look can be
pointed to with :code:`-DOpenMP_ROOT=<folder>`: Homebrew's :code:`libomp` on
macOS needs this, e.g. :code:`-DOpenMP_ROOT=$(brew --prefix libomp)`. (With
:code:`-DUSE_OPENMP=ON` alone CMake prints a warning instead and builds the
libraries without OpenMP.) The configure output reports what is built in the
lines ``OpenMP parallelization is ON`` and ``Parallel PSD projection is ON``.

For the Python package the option is :code:`use_parallel_psd`, in an install
from source: see :ref:`the Python install page <python_openmp_install>`.

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
thread takes a few microseconds to start. SCS counts 2 microseconds for each
thread, and 5 on macOS (:code:`PSD_PAR_THREAD_COST` in :code:`src/cones.c`).
The number of threads grows with about the square root of the time that SCS
estimates for the projection on one thread. With many blocks of similar size,
a projection of less than about 13 microseconds (32 on macOS) stays on one
thread, one of 50 microseconds gets about 6 threads (4 on macOS), and one of a
millisecond gets up to about 30 (20 on macOS). A solve uses fewer threads if
the OpenMP thread count has been lowered since the setup.

Even with the option, the blocks are projected one after another, and the BLAS
is free to use its own threads, when

* the library was compiled without OpenMP (:code:`scsmkl` and :code:`scscudss`
  built with CMake),
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

The number of threads is chosen as if the threads were awake when a projection
starts. With the GNU OpenMP runtime on Linux they wait actively for a short
time after a parallel region (:code:`GOMP_SPINCOUNT`) and then sleep, and
waking them can cost more than a small projection gains. On one machine, PSD
blocks that took less than a millisecond in all made an iteration up to 27%
slower when the rest of the iteration took a millisecond or more. If the PSD
blocks are a small part of a large problem, a build without
:code:`USE_PARALLEL_PSD` suits it better. :code:`OMP_WAIT_POLICY=active` keeps
the threads awake between projections, at the price of cores that are busy
while they wait.

What SCS checks, and where it cannot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCS guards against two of the cases in the table, on some platforms only. That
is the reason the option is left to whoever builds SCS and knows its BLAS.

* MKL on another OpenMP runtime: SCS sets the thread-local thread count of MKL
  to one during the projection. This goes through a weak reference to
  :code:`MKL_Set_Num_Threads_Local`, which needs an ELF platform such as Linux
  and the MKL symbols visible to SCS. It does not work on Windows or macOS,
  with a static MKL linked without the MKL Pardiso linear solver (which
  references the function and so links it in), or with an MKL loaded
  privately by another library. A static MKL on its GNU OpenMP layer in a
  build with :code:`USE_OPENMP` and GCC (the scs-python option
  :code:`mkl_static_threading=gnu_thread`, see
  :ref:`python_mkl_threaded_install`) does not depend on the guard: SCS and MKL
  share one runtime, which runs MKL on one thread inside SCS's parallel
  region.
* A pthreads build of OpenBLAS with more than one thread: on ELF platforms SCS
  checks for it at setup and then projects the blocks one after another. With
  :code:`OPENBLAS_NUM_THREADS=1` the projection is parallel, but the BLAS calls
  of the threads still wait for one lock inside an OpenBLAS that was not built
  with OpenMP. SCS counts that when it chooses the number of threads: small
  blocks mostly stay serial, and blocks of :math:`10 \times 10` to
  :math:`32 \times 32` get two to six threads. An OpenMP build of OpenBLAS
  does not have this limit. Both checks need the functions
  :code:`openblas_get_parallel` and :code:`openblas_get_num_threads` visible to
  SCS: they are off for a static :code:`libopenblas.a` (the weak references
  can stay unresolved), for an OpenBLAS with suffixed symbol names, and for
  one loaded privately by another library (libblastrampoline in Julia).
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
