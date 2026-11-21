.. _blas_lapack:

BLAS and LAPACK
===============
BLAS and LAPACK are dense numerical linear algebra packages. SCS uses these
libraries in the following places:

* To compute the eigen-decomposition for the semidefinite :ref:`cone <cones>`
  projection (see :ref:`below <blas_lapack_psd>`), and the
  eigen-decompositions and singular value decompositions for the
  :ref:`spectral cones <spectral_cones>`
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


.. _blas_lapack_psd:

Projection onto the PSD cone
----------------------------
SCS projects each positive semidefinite (PSD) block onto the cone with a LAPACK
eigendecomposition, and rebuilds the block from the eigenvectors of its
positive eigenvalues with a BLAS product (:code:`syrk`, or :code:`herk` for a
complex block). For a real block the eigendecomposition is LAPACK's
divide-and-conquer routine :code:`syevd` (for a complex block see
:ref:`blas_lapack_psd_complex`). On most SDPs this is where most of the time
of an iteration goes, so the LAPACK library, and the
:ref:`number of threads <blas_lapack_threads>` it runs on, matter.

What it costs and gains
~~~~~~~~~~~~~~~~~~~~~~~
Earlier versions of SCS used the MRRR routine :code:`syevr`. Compared with it,
:code:`syevd`:

* **is faster**: on the SDPLIB problems the projection takes 29% to 51% less
  time per iteration (geometric mean over the problems) with OpenBLAS, NVPL,
  MKL and Apple Accelerate, and a whole iteration 28% to 47% less (OpenBLAS,
  NVPL and Accelerate). It gains most on blocks with exactly repeated
  eigenvalues, which :code:`syevr` is slow on: ``thetaG11`` is projected 2.7
  to 2.9 times faster. With the BLAS on one or four threads, no problem's
  projection took more than 2.3% longer per iteration, except ``qap8`` on a
  slow core of a chip with fast and slow cores (11% to 16%).
* **is as accurate or more**: errors of :math:`2` to :math:`4 \times 10^{-15}`
  relative to the largest entry of the projection, against :math:`3` to
  :math:`9 \times 10^{-15}`, on blocks with repeated eigenvalues.
* **needs more memory**: :code:`syevd`'s work array has about :math:`2n^2`
  entries for a block of size :math:`n`, where :code:`syevr` needed
  :math:`n^2` for the eigenvectors, so the projection of a block takes about
  :math:`3n^2` floating-point numbers instead of :math:`2n^2`. In double
  precision that is 600 MB instead of 400 MB for a block of size 5000.
* **starts the threads of a multithreaded BLAS at smaller blocks**: with the
  BLAS at its default, on all 20 cores of a chip with 10 fast and 10 slow
  cores, problems whose largest block has size 45 to 130 took up to 2.1 times
  as long in the projection as with :code:`syevr` (``control7``), and
  OpenBLAS kept 15 to 20 cores busy where :code:`syevr` kept 3 to 12. Such
  problems are faster on one thread than at the default, with either routine:
  see :ref:`blas_lapack_threads`.

.. _blas_lapack_psd_complex:

Complex blocks
~~~~~~~~~~~~~~
SCS decomposes a complex (Hermitian) PSD block with :code:`heevd`, the complex
counterpart of :code:`syevd`, where earlier versions used :code:`heevr`. We
compared the two in whole solves of 28 problems with complex blocks of size 2
to 300 (quantum fidelity, PPT relaxations, channel capacity, PhaseLift, phase
synchronization, complex max-cut and SCS's tests), with Apple Accelerate, MKL,
NVPL and OpenBLAS. Both solved every problem, in the same number of
iterations, and over all of them took the same time (:code:`heevr` 0.3% less,
shifted geometric mean). Compared with :code:`heevr`, :code:`heevd`:

* **is faster with MKL, NVPL and OpenBLAS**: :code:`heevr` took 1.15 to 1.25
  times as long in the projection per iteration, and with MKL longer at every
  block size.
* **is faster on blocks with nearly repeated eigenvalues**: on the PPT
  relaxations, where about half of the eigenvalues come in pairs closer than
  :math:`10^{-8}` relative to each other, the solves took 1.07 to 1.69 times
  as long with :code:`heevr`, with every library.
* **is faster on small blocks**: :code:`heevr` took 1.07 to 1.43 times as
  long in the projection on blocks of size 64 or less, with every library.
* **is slower on larger blocks with Accelerate**: there :code:`heevr` took
  0.91 of the projection time on blocks of size 65 to 150 and 0.71 to 0.74
  above. Problems with a block of size 200 to 300 took 1.2 to 1.5 times as
  long to solve with :code:`heevd`, and all the problems together 1.1 times
  (shifted geometric mean). With NVPL, :code:`heevr` took 0.91 of the
  projection time on blocks larger than 150.
* **needs more memory**: the projection of a complex block of size :math:`n`
  takes about :math:`6n^2` floating-point numbers instead of :math:`4n^2`.

SDPLIB has no complex blocks.

If LAPACK does not converge
~~~~~~~~~~~~~~~~~~~~~~~~~~~
When LAPACK reports that the eigendecomposition of a block did not converge,
SCS decomposes the block again with the QR-based :code:`syev` (:code:`heev`
for a complex block). That is slower: forced on 11 SDPLIB problems, the
projection took 1.4 to 3.6 times as long, and the problems were solved as with
:code:`syevd`. It has happened with Apple Accelerate, whose :code:`dsyevd` did
not converge on one iterate of ``thetaG11`` (a block of size 801), and whose
:code:`dsyevr` did not converge on the block of
`issue #531 <https://github.com/cvxgrp/scs/issues/531>`_. If the QR-based
routine fails as well, the solve stops with status ``failure`` and the message
``error in project_cones``.

Very large blocks
~~~~~~~~~~~~~~~~~
LAPACK counts the entries of its work arrays in its own integers. With 32-bit
integers, the usual case, the :math:`1 + 6n + 2n^2` entries that
:code:`syevd` needs fit up to :math:`n = 32766`. For a larger block, the setup
fails with a message saying so. Use a BLAS and LAPACK with 64-bit integers
(the :ref:`compiler flag <compile_flags>` :code:`BLAS64`) for blocks that
large, and for blocks larger than :math:`n = 46340` SCS's own 64-bit
integers as well (the compiler flag :code:`DLONG`): SCS indexes the entries
of a block with numbers up to :math:`n^2`, which then do not fit in 32 bits.


.. _blas_lapack_threads:

Number of BLAS threads
----------------------
The :ref:`projection <blas_lapack_psd>` of a PSD block is a LAPACK
eigendecomposition (:code:`syevd` for a real block) and a BLAS product
(:code:`syrk` or :code:`herk`). A multithreaded BLAS can split these calls,
and SCS's other BLAS calls (see the list above), over its own threads, by
default about one per core. Waking and synchronizing the threads has a fixed
cost. In our measurements the default made the projection of blocks of size
up to 300 slower than one thread did, with both libraries; a few threads
gained on blocks of size 500 and more, and on smaller blocks it depended on
the library. The number of threads is set by an environment variable of the
BLAS, so you can choose.

When to set one thread
~~~~~~~~~~~~~~~~~~~~~~
Put the BLAS on one thread for an SDP whose PSD blocks are all of size 300 or
less, in particular when several solves run at the same time, since each would
start threads for every core. At their default, the two libraries we measured
took 1.5 to 2.5 times as long as on one thread on blocks of size 37 to 300
(geometric means, see :ref:`blas_lapack_threads_measured`). On 4 threads NVPL
was slower than on one as well (1.18 to 1.28 times), while OpenBLAS was faster
on blocks of size 37 and more (0.69 to 0.88 of the time): with OpenBLAS,
compare one thread with a few.

Use a few threads instead when

* the problem has PSD blocks of size 500 or more. On 4 fast cores the
  projection took 0.85 (NVPL) and 0.47 (OpenBLAS) of its time on one thread at
  size 500, and 0.58 and 0.41 at sizes 800 to 2000. The default, all 20 cores
  of that machine, was slower than one thread at size 500 with NVPL (1.52),
  and with OpenBLAS gained less than 4 threads did (0.89 against 0.47), so
  compare a few thread counts, and on a machine with fast and slow cores count
  only the fast ones;
* the linear system is large and solved with MKL Pardiso, which runs on MKL's
  threads (see the threads paragraph of :ref:`MKL Pardiso <mkl>`), or with the
  :ref:`dense direct <dense>` solver, which forms and factorizes its matrix
  with BLAS and LAPACK (not measured).

A problem with several PSD blocks can instead use threads across the blocks,
with the BLAS on one thread inside each, in a build with the parallel PSD
projection (off by default): see :ref:`blas_lapack_openmp`.

Which variable to set
~~~~~~~~~~~~~~~~~~~~~
Set it to 1, or to the number of threads to use.

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Your BLAS
     - Set
     - Notes
   * - OpenBLAS with a thread pool of its own: built with pthreads, its
       default (the default :code:`libopenblas-dev` of Debian and Ubuntu is
       one), or on Windows with Windows threads
     - :code:`OPENBLAS_NUM_THREADS=1`
     - When that is unset it reads :code:`OMP_NUM_THREADS`, so
       :code:`OMP_NUM_THREADS=1` works too
   * - A BLAS built with OpenMP: an OpenMP build of OpenBLAS, NVPL's OpenMP
       (:code:`_gomp`) libraries
     - :code:`OMP_NUM_THREADS=1`
     - An OpenMP build of OpenBLAS ignores :code:`OPENBLAS_NUM_THREADS`
   * - MKL
     - :code:`MKL_NUM_THREADS=1`
     - Takes precedence over :code:`OMP_NUM_THREADS`. It also puts MKL Pardiso
       on one thread. A sequential MKL, such as the MKL backend of the x86-64
       Linux Python wheels (their default linear solver), always runs on one
       thread. In the projection, MKL 2026.1 on 4 threads of an AMD Zen 3
       took 1.8 to 7 times as long as on one thread for real blocks of size 2
       to 10 and 1.1 to 1.9 times for size 16 to 200; from size 256 on its
       threads were faster
   * - Accelerate (macOS)
     - Nothing
     - :code:`VECLIB_MAXIMUM_THREADS=1` did not change the time of the
       projection on an Apple M4 (82 SDPLIB problems: time at the default over
       time with the variable 0.96 to 1.07, geometric mean 1.00)

Set the variable before the program starts (in Python, at the latest before
importing NumPy or SCS): the libraries read it when they start. It holds for
every call to that library in the process, including SCS's dense direct
linear solver and your own code.

In a build of SCS with OpenMP (:code:`USE_OPENMP`; the pre-built Python wheels
are built without it), :code:`OMP_NUM_THREADS` is also SCS's own thread count,
and setting it to 1 stops SCS's own threads too, the :ref:`parallel PSD
projection <blas_lapack_openmp>` included. With OpenBLAS's own thread pool or
with MKL, set the library's variable instead. An OpenMP build of OpenBLAS, or
NVPL, has no variable of its own, so :code:`OMP_NUM_THREADS` is one budget for
SCS and the BLAS. Inside a parallel PSD projection the BLAS runs on one thread
by itself. Where a build with :code:`USE_PARALLEL_PSD` projects the blocks one
after another (one block that holds most of the work, say), it sets the OpenMP
thread count, and MKL's thread-local count, to one for the blocks of size 384
or less, so such a BLAS and MKL run on one thread there too, and on all the
threads only for larger blocks (see :ref:`the serial projection
<blas_lapack_openmp_serial>`). In a build without
:code:`USE_PARALLEL_PSD` the BLAS runs on all the threads for every block. We
measured that case in a build of SCS without OpenMP: on control7 (blocks of
size 70 and 35), NVPL with :code:`OMP_NUM_THREADS=4` made the projection take
1.82 times as long as with 1, and a whole iteration 1.37 times. For such a
problem compare :code:`OMP_NUM_THREADS=1` with your thread count.

To check, watch the CPU use of the process during a solve (``top``). With the
BLAS on one thread and SCS built without OpenMP it stays near 100%; several
hundred percent on a problem with small PSD blocks means that the BLAS threads
are running.

.. _blas_lapack_threads_measured:

What we measured
~~~~~~~~~~~~~~~~
Time of the PSD projection (with :code:`syevd`) per iteration on the 58
SDPLIB problems whose largest PSD block is of size 37 or more, relative to the
BLAS on one thread (geometric mean, range in brackets; below 1 is faster). The
machine is an NVIDIA GB10 with 10 fast and 10 slow ARM cores, the libraries
NVPL and a pthreads build of OpenBLAS 0.3.31. The runs on one and on 4
threads were pinned to fast cores; 20 threads is the default of both libraries
there.

.. list-table::
   :header-rows: 1
   :widths: 24 19 19 19 19

   * - Size of the largest PSD block
     - NVPL, 4 threads
     - NVPL, 20 threads
     - OpenBLAS, 4 threads
     - OpenBLAS, 20 threads
   * - 37 to 124 (25 problems)
     - 1.28 [1.04, 2.12]
     - 2.55 [1.26, 9.43]
     - 0.88 [0.81, 0.96]
     - 1.50 [1.07, 2.24]
   * - 150 to 300 (17)
     - 1.18 [1.06, 1.33]
     - 2.16 [1.98, 2.46]
     - 0.69 [0.61, 0.79]
     - 1.52 [0.93, 2.07]
   * - 500 (8)
     - 0.85 [0.84, 0.86]
     - 1.52 [1.48, 1.58]
     - 0.47 [0.45, 0.49]
     - 0.89 [0.66, 0.96]
   * - 800 to 2000 (8)
     - 0.58 [0.50, 0.66]
     - 0.79 [0.47, 1.14]
     - 0.41 [0.36, 0.46]
     - 0.45 [0.29, 0.69]

* control7 (blocks of size 70 and 35) took 0.39 ms per projection with NVPL
  on one thread, 0.70 ms on 4 and 3.65 ms on 20, and a whole iteration
  1.14 ms, 1.56 ms and 5.10 ms.
* The pthreads OpenBLAS on 4 threads was faster than on one on every problem
  with blocks of size 37 or more. At its default it was slower up to size 300
  (1.50 times at sizes 37 to 124, 1.52 at 150 to 300) and kept cores busy for
  nothing: 9 problems with blocks of size 37 to 124 used 15 to 20 cores
  (1526% to 1998% CPU) where one thread was 1.07 to 2.17 times as fast, and
  truss8, whose blocks are of size 19 at most, used 337%. Those cores are
  taken from anything else on the machine, such as other solves.
* The 20-thread columns include the slow cores and depend on where the system
  placed the threads: the 27 problems whose blocks are all smaller than 37
  were up to 2.4 times slower there, though 4 threads changed their time by
  at most 10%. So they show what the default can cost on such a machine, not
  the cost of the threads alone.
* We did not measure x86 (MKL, or OpenBLAS on x86) or a machine with only one
  kind of core, where the sizes at which threads start to pay may differ.
  There, compare the time of a solve on one and on a few threads.


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
       times faster than the serial one on 4 threads (measured before SCS moved
       to :code:`syevd`)

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

SCS estimates the time of each block from its size, with a formula
(:code:`psd_block_cost` in :code:`src/cones.c`) fitted to the time
:code:`syevd` takes for a real block with NVPL on one core of an Arm
Cortex-X925. The estimate decides how the blocks are shared out and, with
the time per thread, how many threads are used. BLAS libraries differ most on
small blocks: Accelerate on an Apple M4 and OpenBLAS take a quarter to two
thirds of NVPL's time for blocks of size 3 to 10, and MKL on x86-64 (AMD
Zen 3) 1.1 times. SCS uses one formula for all of them. Following NVPL, the
slowest of the first three, was the fastest of the formulas we compared over
those three libraries, but not with each of them. On the block structures
whose split changed, against the previous formula, the projection took:

* with NVPL, 0.84 times the time at 2 to 4 threads and 0.89 at 10 threads;
* with OpenBLAS, 0.94 and 1.20 times (see the OpenBLAS lock below);
* with Accelerate, 1.06 and 1.12 times;
* with MKL (static, on the GNU OpenMP runtime, as the scs-python option
  :code:`mkl_static_threading=gnu_thread` builds it), measured afterwards:
  0.75, 0.79 and 0.73 times at 2, 3 and 4 threads (10 threads not measured).
  Part of that gain: in this build, blocks projected one after another each
  run on MKL's own threads, which took longer than one thread at every block
  size measured (:math:`3 \times 3` to :math:`128 \times 128`).

Whole solves of the SDPs among those structures reached the same results in
0.98 times the time overall (0.86 with NVPL, 1.06 with OpenBLAS, 1.07 with
Accelerate), and 1.08 times at 10 threads alone. With MKL they reached the
same verified results in 0.81, 0.59 and 0.62 times the time at 2, 3 and 4
threads. The cost is on the libraries that are faster on small blocks: a
projection of 9 to 62 microseconds made of 5 to 64 blocks of size 3 to 10 can
get threads, or more of them, where fewer are faster (up to 2.3 times on the
M4). If that matters, compare the time of a solve with a lower
:code:`OMP_NUM_THREADS`. OpenBLAS on x86-64 was not measured.

Complex blocks keep the estimate fitted to :code:`heevr`, which is close to
:code:`heevd`'s time with NVPL up to size 100 and up to 1.5 times too low for
larger blocks: the estimate fitted to :code:`heevd` was slower in the same
measurements.

Even with the option, the blocks are projected one after another, and the BLAS
is free to use its own threads (for blocks larger than :math:`384 \times 384`,
see :ref:`below <blas_lapack_openmp_serial>`), when

* the library was compiled without OpenMP (:code:`scsmkl` and :code:`scscudss`
  built with CMake),
* the threads would gain little: the blocks are too little work, or one block
  holds most of it (the parallel projection is no faster than its largest
  block),
* a block larger than :math:`512 \times 512` holds more than half of the work
  (the parallel projection would then gain less than a factor of two, and a
  multithreaded BLAS may do better on that block). Whether it does depends on
  the BLAS. MKL's threads beat the parallel projection on such a block larger
  than 384 (AMD Zen 3, 3 and 4 threads), NVPL's from about 1000 (DGX Spark),
  and Accelerate's never (Apple M4). So where this rule projects such a block
  of size 640 to 2000 one after another, it took 1.03 to 1.85 times the time
  of the parallel projection on Accelerate, and 1.02 to 1.56 times with a BLAS
  on one thread (a pthreads OpenBLAS with :code:`OPENBLAS_NUM_THREADS=1`).
  512 lies between the libraries: over all four we measured (with OpenBLAS's
  OpenMP build) it took 1.6% more time than the best size, 768, with the
  smaller blocks held to one thread, and 2.4% more than the best, 1024 and
  1500, with them left to the BLAS's threads,
* the BLAS is a pthreads build of OpenBLAS set to more than one thread, where
  SCS can see that (see below), or
* the solver is set up inside an active OpenMP parallel region, where the
  threads of the caller are the parallelism.

.. _blas_lapack_openmp_serial:

When a build with the option projects the blocks one after another, it sets
the OpenMP thread count to one for each block of size :math:`2 \times 2` to
:math:`384 \times 384`, and back to what it was after them
(:code:`PSD_SERIAL_HOLD_N` in :code:`src/cones.c`). Larger blocks run with the
thread count as it was. This reaches a BLAS on the OpenMP runtime of SCS that
takes its number of threads from the OpenMP thread count, such as NVPL and an
OpenMP build of OpenBLAS. For MKL, whose own thread count takes precedence
over the OpenMP one, it also sets MKL's thread-local count to one for the same
blocks, as the parallel projection does. That goes through a weak reference to
MKL's symbol, which needs Linux (ELF) and MKL's symbols visible to SCS. Such a
BLAS, and MKL, then run these blocks on one thread, as inside a parallel
projection, because on them their threads cost more than they gain:

* MKL 2026.1 (statically linked, GNU OpenMP layer, AMD Zen 3: the CPU build of
  the SCS benchmark, which sets :code:`MKL_NUM_THREADS`): blocks of size 128 or
  less took the time on one thread, where MKL's 4 threads took 1.8 to 7 times
  as long for size 2 to 10. Whole solves of the 148 pairs of problem and
  tolerance among the benchmark's SDPs whose blocks are projected one after
  another, at 4 threads: 127 pairs solved against 126, in 0.95 of the wall
  clock (shifted geometric mean, shift 10 s; 0.85 at 10 ms). MKL's threaded
  eigensolver rounds differently from one thread, so the iterates change, and
  the pairs that changed outcome went both ways (5 solved only with the hold,
  4 only without it).

* NVPL (DGX Spark): real blocks of size 32 to 300 took 0.78 to 0.95 of the
  time at 2 and 4 threads, and real blocks of size 32 to 384 0.48 to 0.87 of
  it at 10 threads. Whole solves of the 70 shorter pairs of problem and
  tolerance among the benchmark's SDPs whose blocks are projected one after
  another (66 pairs, and 4 projected in parallel, where nothing is held)
  reached the same results, with the same iterates, in 0.94 of the wall clock
  at 4 threads and 0.82 at 10 (shifted geometric mean, shift 10 s; 0.93 and
  0.78 at 10 ms).
* An OpenMP build of OpenBLAS: 0.59 to 0.98 of the time for blocks of size 32
  to 384 at 4 threads on an Apple M4, and 0.15 to 0.82 for real blocks of size
  3 to 512 at 4 threads on an AMD Zen 3.

Where the OpenMP thread count does not reach the BLAS, it changes nothing,
neither the results nor the time:

* a BLAS with threads of its own: a pthreads build of OpenBLAS (the default on
  Debian and Ubuntu, the SciPy and NumPy wheels, Julia) and Accelerate;
* a BLAS on another OpenMP runtime than that of SCS;
* an OpenMP build of OpenBLAS from version 0.3.34 on, once the program has
  called :code:`openblas_set_num_threads`: it then keeps that count (from its
  source; not measured);
* MKL where SCS cannot see it through the weak reference (macOS, Windows, or
  an MKL loaded privately by another library, as Julia's libblastrampoline
  does): its own thread count takes precedence over the OpenMP one. MKL's
  threads are slow on small blocks (see the MKL row of the table in
  :ref:`blas_lapack_threads`). If your problem's blocks are small and
  projected one after another, :code:`MKL_NUM_THREADS=1` puts such an MKL on
  one thread, but for the whole process, MKL Pardiso included;
* a sequential BLAS, and builds of SCS without OpenMP (the Python wheels).

We did not measure MKL without :code:`MKL_NUM_THREADS`, MKL on Intel's OpenMP
runtime, BLIS, or Windows. The cost is on the larger blocks below the size:
there the threads of the BLAS can pay, and the block runs on one thread. A
block of size 512 took 1.10 to 1.20 times as long on one thread as on 4 (NVPL,
OpenBLAS on the M4), and blocks of size 200 to 384 1.02 to 1.12 times as long
as on 2 threads (OpenBLAS on the M4). Where the hold works, sizes 384 and 512
came out the same over all the block sizes and thread counts we measured, and
384 risks less. With LLVM's OpenMP runtime (libomp, as clang links it),
setting the thread count to one also shrinks the team of the calling thread,
which the next parallel region grows again: about 1 microsecond per
projection at 4 threads on an Apple M4, in a microbenchmark. GNU's runtime
does not show this.

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
  with OpenMP. SCS counts that when it chooses the number of threads: a few
  small blocks stay serial, many blocks of :math:`3 \times 3` to
  :math:`10 \times 10` get two to four threads, and blocks of
  :math:`16 \times 16` to :math:`32 \times 32` four to six. It does not
  count the cores the threads run on. On a DGX Spark (two clusters of Arm
  Cortex-X925 cores), a team of the same size took about 1.47 times as long
  with 10 OpenMP threads on the fast cores of both clusters as with 4 on one
  cluster (geometric mean of 55 ratios, one per block structure and team
  size, over 44 structures; 1.31 times over 11 ratios and 10 structures with
  NVPL, which has no lock). The two thread counts were timed in separate
  sessions, not side by side. With 10 threads the projection took 1.20 times
  the time of the previous formula on the block structures whose split
  changed, against 0.94 times at 2 to 4 threads. An OpenMP build of OpenBLAS
  does not have this limit. Both checks need the functions
  :code:`openblas_get_parallel` and :code:`openblas_get_num_threads` visible
  to SCS: they are off for a static :code:`libopenblas.a` (the weak
  references can stay unresolved), for an OpenBLAS with suffixed symbol
  names, and for one loaded privately by another library (libblastrampoline
  in Julia).
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

A parallel projection takes more memory: each thread has its own
eigendecomposition workspace, sized for the blocks that thread projects, and
none is larger than the workspace of the serial projection. The extra memory
is therefore at most the PSD workspace of the serial projection times the
number of threads less one, and less when only a few threads get the largest
blocks. That workspace is about :math:`3n^2` floats for the largest real
block of size :math:`n` (see :ref:`above <blas_lapack_psd>`, and
:ref:`blas_lapack_psd_complex` for a complex block).

Each block goes through the same LAPACK calls, with the same work-array sizes,
as in the serial projection, so with single-threaded BLAS the result is the
same bit for bit, whatever the number of threads (unless the rounding of the
BLAS depends on the memory alignment of its arguments, as it can with MKL).
