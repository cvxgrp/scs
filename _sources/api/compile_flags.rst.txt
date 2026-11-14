.. _compile_flags:

Compile flags
-------------

Typically the user will not need to interact with these flags, but we document them
here for completeness.  The flags in the table below are defined in
:code:`scs.mk` and can be overridden via the command line when compiling with
the Makefile by executing, e.g., :code:`make DLONG=1`, to set the :code:`DLONG`
flag to True.  Only some of them are available in the CMake build, see
:ref:`below <compile_flags_cmake>`.  Note this list is not exhaustive: more
variables of the Makefile are listed :ref:`below <compile_flags_make>` (those
of the CUDA builds are on the :ref:`install <c_install>` page), and further
constants are not exposed but defined in :code:`include/glbopts.h`.


.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Name
     - Description
     - Permitted values
     - Default
   * - :code:`DLONG`
     - If True use 64 bit integers, else 32 bit
     - True/False
     - 0
   * - :code:`SFLOAT`
     - If True use 32 bit floats, else 64 bit. WARNING: experimental. Single precision is not covered by the test suite, and most of the solver tests fail when built with it. It does not work with Apple's Accelerate BLAS, which is what :code:`make` links on macOS when :code:`USE_LAPACK=1`: Accelerate returns a :code:`double` from the single precision :code:`dot` and :code:`nrm2`, where SCS expects a :code:`float` (see :code:`src/linalg.c`)
     - True/False
     - 0
   * - :code:`CTRLC`
     - Listen to CTRL-C interruptions
     - True/False
     - 1
   * - :code:`NO_TIMER`
     - Disables code timing
     - True/False
     - 0
   * - :code:`NO_VALIDATE`
     - Disables data validation
     - True/False
     - 0
   * - :code:`NO_PRINTING`
     - Disables all printing in compiled binary
     - True/False
     - 0
   * - :code:`NO_READ_WRITE`
     - Disables the read/write code
     - True/False
     - 0
   * - :code:`GPU_TRANSPOSE_MAT`
     - If on GPU store A transpose in memory
     - True/False
     - 1
   * - :code:`VERBOSITY`
     - Verbosity level (for debugging)
     - :math:`\mathbf{N}`
     - 0
   * - :code:`USE_LAPACK`
     - Whether to link in :ref:`BLAS/LAPACK <blas_lapack>`
     - True/False
     - 1
   * - :code:`USE_OPENMP`
     - Use openmp to parallelize some computation (matrix-vector products and the projection onto exponential cones). The OpenMP thread count (:code:`OMP_NUM_THREADS`) is the number of threads
     - True/False
     - 0
   * - :code:`USE_PARALLEL_PSD`
     - Project the :ref:`PSD blocks <blas_lapack_openmp>` of a problem in parallel. Requires :code:`USE_OPENMP=1` and :code:`USE_LAPACK=1`. Only faster when the BLAS runs on one thread inside the threads of SCS, so it is off by default
     - True/False
     - 0
   * - :code:`USE_SPECTRAL_CONES`
     - Enable spectral cone support (log-determinant, nuclear norm, :math:`\ell_1` norm, sum-of-largest-eigenvalues). Requires :code:`USE_LAPACK=1`.
     - True/False
     - 0
   * - :code:`BLAS64`
     - BLAS/LAPACK uses 64-bit integers (ILP64)
     - True/False
     - 0
   * - :code:`NOBLASSUFFIX`
     - The BLAS library has no function name suffix
     - True/False
     - 0
   * - :code:`BLASSUFFIX`
     - The BLAS library uses this suffix
     - String
     - '_'
   * - :code:`MATLAB_MEX_FILE`
     - If compiling for use in MATLAB
     - True/False
     - 0
   * - :code:`PYTHON`
     - If compiling for use in python
     - True/False
     - 0
   * - :code:`USING_R`
     - If compiling for use in R: SCS prints with :code:`Rprintf`, so the R headers must be on the include path. It defines the macro :code:`R_LANG`, which is the name the source reads and what the R package passes in :code:`CFLAGS`. Build only the libraries, e.g., :code:`make out/libscsdir.a USING_R=1 CFLAGS=-I<R include folder>`: the demo and test programs call :code:`Rprintf` and do not link without R
     - True/False
     - 0

.. _compile_flags_make:

Other Makefile variables
^^^^^^^^^^^^^^^^^^^^^^^^

The Makefile also uses the variables below. They are set in the same way, e.g.,
:code:`make OPT=-O3 BLASLDFLAGS=-lopenblas`.

.. list-table::
   :widths: 20 55 25
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - :code:`OPT`
     - Optimization flags for the compiler. The default includes :code:`-march=native`, so the libraries are tuned to the processor of the machine they are built on and may not run on a different one. Use, e.g., :code:`OPT=-O3` for libraries that will run on other machines. (The CMake build does not use :code:`-march=native`)
     - :code:`-O3 -march=native -fno-math-errno`
   * - :code:`BLASLDFLAGS`
     - Linker flags for the :ref:`BLAS and LAPACK <blas_lapack>` libraries, e.g., :code:`BLASLDFLAGS=-lopenblas` or :code:`BLASLDFLAGS="-L<folder> -lopenblas"`
     - :code:`-llapack -lblas` (empty if :code:`USE_LAPACK=0`)
   * - :code:`PREFIX`
     - Where :code:`make install` copies the header files (to :code:`$(PREFIX)/include/scs`) and the direct and indirect libraries (to :code:`$(PREFIX)/lib`). :code:`DESTDIR`, if set, is put in front of :code:`PREFIX`
     - :code:`/usr/local`
   * - :code:`MKLROOT`
     - Location of MKL, for :code:`make mkl`
     - Not set
   * - :code:`OMPROOT`
     - Folder whose :code:`lib` subfolder holds the Intel OpenMP library (:code:`libiomp5`) that :code:`make mkl` links, typically :code:`/opt/intel/oneapi/compiler/latest`
     - :code:`$(MKLROOT)`
   * - :code:`COVERAGE`
     - If True compile with :code:`--coverage` to generate test coverage data
     - 0
   * - :code:`CFLAGS`
     - Additional compiler flags, e.g., :code:`CFLAGS=-DSPECTRAL_DEBUG`. SCS adds its own flags after them
     - Empty

:code:`OPT`, :code:`PREFIX`, :code:`MKLROOT`, :code:`OMPROOT` and :code:`CFLAGS`
can also be set as environment variables. :code:`BLASLDFLAGS`, :code:`COVERAGE`
and the flags of the first table are only read from the command line of
:code:`make`.

.. _compile_flags_cmake:

CMake
^^^^^

With CMake the flags are options that are set when configuring, e.g.,
:code:`cmake -DDLONG=ON ..`. The following flags of the first table are CMake
options, with the same default: :code:`DLONG`, :code:`SFLOAT`,
:code:`NO_PRINTING`, :code:`NO_READ_WRITE`, :code:`USE_LAPACK`, :code:`BLAS64`,
:code:`USE_OPENMP`, :code:`USE_PARALLEL_PSD` and :code:`USE_SPECTRAL_CONES`.
:code:`GPU_TRANSPOSE_MAT` is a CMake option too, but only together with
:code:`-DUSE_GPU=ON`.

:code:`CTRLC` is always enabled in a CMake build. :code:`NO_TIMER`,
:code:`NO_VALIDATE`, :code:`VERBOSITY`, :code:`NOBLASSUFFIX`,
:code:`BLASSUFFIX`, :code:`MATLAB_MEX_FILE`, :code:`PYTHON` and :code:`USING_R`
are not CMake options, and CMake ignores them with a warning if they are passed
as one (:code:`-DNO_TIMER=ON`, say). They can be passed to the compiler
directly, e.g., :code:`cmake -DCMAKE_C_FLAGS="-DNO_TIMER=1" ..` (for R the
macro is :code:`R_LANG`, so :code:`-DCMAKE_C_FLAGS="-DR_LANG"`).

CMake finds the BLAS and LAPACK libraries itself. To choose a particular
library set CMake's `BLA_VENDOR
<https://cmake.org/cmake/help/latest/module/FindBLAS.html>`_, e.g.,
:code:`cmake -DBLA_VENDOR=OpenBLAS ..`. If the :code:`MKLROOT` environment
variable is set and :code:`BLA_VENDOR` is not, SCS asks CMake for MKL. If CMake
finds no LAPACK, SCS links :code:`blas` and :code:`lapack` by name, except with
:code:`BLAS64=ON`, where the configuration stops with an error.

If CMake finds no OpenMP, :code:`USE_OPENMP=ON` gives a warning and libraries
without OpenMP (the configure output then says ``OpenMP parallelization is
OFF``), except with :code:`USE_PARALLEL_PSD=ON`, where the configuration stops
with an error. OpenMP that is installed where CMake does not look can be
pointed to with :code:`-DOpenMP_ROOT=<folder>`, which Homebrew's
:code:`libomp` on macOS needs.

The options that select what is built (:code:`BUILD_SHARED_LIBS`,
:code:`BUILD_TESTING` and the options for the other linear solvers) are
described in the :ref:`install <c_install>` and
:ref:`linear solver <linear_solver>` pages. CMake has two more options, which
the Makefile does not have (there, pass the definition in :code:`CFLAGS`):

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Name
     - Description
     - Permitted values
     - Default
   * - :code:`USE_SPECTRAL_TIMING_FLAG`
     - Time the projections onto the matrix cones and the spectral vector cones. Requires :code:`USE_SPECTRAL_CONES=ON`. Defines :code:`SPECTRAL_TIMING_FLAG`, which adds the fields :code:`ave_time_matrix_cone_proj` and :code:`ave_time_vector_cone_proj` (average time per iteration, in milliseconds) to the end of :ref:`ScsInfo <ScsInfo>`. Code that calls a library built this way must be compiled with :code:`-DSPECTRAL_TIMING_FLAG` as well, so that both use the same :code:`ScsInfo`
     - ON/OFF
     - OFF
   * - :code:`USE_SPECTRAL_DEBUG`
     - Print debugging output from the spectral cone projections (defines :code:`SPECTRAL_DEBUG`)
     - ON/OFF
     - OFF
