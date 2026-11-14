.. _python_install:

Python
======

The easiest way to install the python version is using `pip <https://pypi.org/project/pip/>`_:

.. code:: bash

  pip install scs

.. important::

   On x86-64 Linux the pre-built manylinux (glibc) wheels include the :ref:`MKL Pardiso <mkl>`
   direct linear solver, linked statically, and SCS selects it
   automatically: it is faster than the built-in solver for most problems,
   often dramatically so on larger ones, and nothing extra needs to be
   installed. Every other wheel falls back to the built-in QDLDL solver.

You can also install directly from source

.. code:: bash

  git clone --recursive https://github.com/bodono/scs-python.git
  cd scs-python
  python -m pip install .

Linear solver backends
----------------------

The pre-built wheels and a from-source install always include two CPU linear
solvers that require no additional dependencies:

- :code:`QDLDL` — the default :ref:`sparse direct solver <direct>` (bundled
  with SCS).
- :code:`CPU_INDIRECT` — the sparse matrix-free solver based on conjugate
  gradients.

They also include :code:`CPU_DENSE`, the :ref:`dense direct solver <dense>`,
unless a source build turns LAPACK off (see below).

The remaining backends require either a platform-specific library (Apple
Accelerate, MKL) or a build-time flag plus an external dependency (CUDA +
cuDSS for GPU). Each section below describes what to install and
how to enable the backend. See :ref:`linear_solver` for an overview of each
solver and :ref:`python_interface` for how to select one at runtime.

Apple Accelerate (macOS)
""""""""""""""""""""""""

On macOS the :ref:`Apple Accelerate backend <apple_accelerate>` is built and
included automatically — no extra install flags are needed. It uses the
Accelerate framework's sparse LDL\ :sup:`T` solver, which is optimized for
Apple hardware including Apple Silicon. The default :code:`linear_solver=scs.LinearSolver.AUTO` selects the
bundled QDLDL on macOS; opt in to Accelerate explicitly with
:code:`linear_solver=scs.LinearSolver.ACCELERATE`.

MKL
"""

On x86-64 Linux the pre-built manylinux wheels include the MKL Pardiso
solver: MKL is linked statically into the ``_scs_mkl`` extension, with the
sequential threading layer and its symbols hidden, so the wheel is
self-contained and the default
:code:`linear_solver=scs.LinearSolver.AUTO` selects MKL with nothing else
installed. MKL is typically faster than QDLDL.

MKL is linked statically rather than bundled as shared libraries, whose
dlopen'd CPU dispatch kernels wheel-repair tools cannot see (see `issue #423
<https://github.com/cvxgrp/scs/issues/423>`_); the backend is
single-threaded, and Intel's license notice ships in the wheel as
``LICENSE-INTEL-MKL.txt``. The aarch64, musllinux, macOS and Windows wheels
and sdist/source builds do not include it.

If your environment already provides MKL (e.g. conda), you can instead build
from source against it:

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_mkl=true .

See :ref:`here <python_interface>` for how to select MKL when solving.

The Windows wheels do not currently include the MKL backend; on Windows MKL
is available in source builds, which use sequential MKL until Intel fixes
the threaded :code:`pkg-config` metadata in its conda packages.

To use 64-bit BLAS/LAPACK integers (ILP64 / :code:`BLAS64`), install with:

.. code:: bash

  python -m pip install -Csetup-args=-Duse_blas64=true -Csetup-args=-Dlink_mkl=true .

The Meson build only supports :code:`BLAS64` together with MKL
(:code:`-Dlink_mkl=true`): ILP64 variants of other BLAS libraries cannot be
verified at build time, and 64-bit prototypes against an LP64 library
corrupt every BLAS call. SCS then requires 64-bit SCS integers as well (the
default in the Meson build). At runtime SCS also checks that the
process-wide MKL interface layer matches the LP64/ILP64 mode it was
compiled for, and fails early if another library already set an
incompatible MKL interface.

Dense direct (LAPACK)
"""""""""""""""""""""

The :ref:`dense direct solver <dense>` reduces the KKT system to a smaller
Gram matrix and factorizes it with LAPACK's Cholesky routines. It is well
suited to small-to-medium problems with a dense constraint matrix :math:`A`,
where dense BLAS/LAPACK outperforms sparse factorization.

It is part of the pre-built wheels and of a default source build: the Meson
option :code:`use_lapack` is on by default, so no extra flag is needed. Only
a source build with :code:`-Csetup-args=-Duse_lapack=false` leaves it out.
Select the backend at runtime with
:code:`linear_solver=scs.LinearSolver.CPU_DENSE`.

Every source build needs a BLAS/LAPACK library that Meson can find, also one
with :code:`use_lapack=false`: Apple Accelerate on macOS, OpenBLAS (or a
reference BLAS and LAPACK) on Linux and Windows, or MKL with
:code:`-Csetup-args=-Dlink_mkl=true`.

.. _python_cudss_install:

GPU direct (cuDSS)
""""""""""""""""""

The :ref:`cuDSS backend <cudss_solver>` runs the sparse direct factorization
and solves on an NVIDIA GPU via `NVIDIA cuDSS
<https://developer.nvidia.com/cudss>`_. On the large LPs and QPs of the
:ref:`benchmarks <benchmarks>` it was faster than SCS on the CPU with MKL
Pardiso. It does not help on SDPs, where the time goes into the cone
projection, which runs on the CPU.

**Prerequisites.** The build links against both the CUDA runtime and cuDSS,
so you need all of the following installed before running
:code:`pip install`:

1. An NVIDIA GPU with a recent CUDA-capable driver.
2. The `CUDA Toolkit <https://developer.nvidia.com/cuda-downloads>`_ (the
   CUDA runtime headers and libraries). Meson finds it itself, from the
   :code:`CUDA_PATH` environment variable (or :code:`CUDA_HOME` or
   :code:`CUDA_ROOT`) or else in :code:`/usr/local/cuda`; it does not look for
   it with :code:`pkg-config`.
3. The `cuDSS library <https://developer.nvidia.com/cudss-downloads>`_. Meson
   finds it with :code:`pkg-config` (:code:`cudss.pc`) or, failing that,
   through the CMake package that cuDSS installs. cuDSS is built for one CUDA
   major version, so use the build that matches your CUDA toolkit: NVIDIA's
   :code:`cudss` apt package installs the builds for CUDA 12 and CUDA 13 and
   makes the CUDA 13 one the default. cuDSS is also available on
   `conda-forge <https://anaconda.org/conda-forge/libcudss>`_ as
   :code:`libcudss` / :code:`libcudss-dev`.

Make sure Meson can find cuDSS: put the directory containing :code:`cudss.pc`
on :code:`PKG_CONFIG_PATH` or, for a cuDSS without :code:`cudss.pc`, set the
environment variable :code:`cudss_DIR` to the folder that holds
:code:`cudss-config.cmake`. Set :code:`CUDA_PATH` if the toolkit is not in
:code:`/usr/local/cuda`, and put the cuDSS shared library on
:code:`LD_LIBRARY_PATH` (Linux) at runtime. For example, with cuDSS unpacked
from NVIDIA's archive under :code:`<cudss-prefix>`:

.. code:: bash

  export cudss_DIR=<cudss-prefix>/lib/cmake/cudss
  export LD_LIBRARY_PATH=/usr/local/cuda/lib64:<cudss-prefix>/lib:$LD_LIBRARY_PATH

**Install.** Once the prerequisites are in place, build SCS with:

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_cudss=true -Csetup-args=-Dint32=true .

The :code:`int32=true` flag is required: the Python build only supports the
GPU solvers with 32-bit integers and stops with an error without it. It
applies to the whole build, so the CPU solvers of that installation use
32-bit integers too. Select the backend at runtime with
:code:`linear_solver=scs.LinearSolver.CUDSS`.

See `this Colab notebook <https://colab.research.google.com/drive/1POCgDNFg8fycHMI9T9N6V3iHFhXRthjn?usp=sharing>`_
for a worked end-to-end example that installs CUDA, cuDSS, and the cuDSS
build of SCS, then solves a problem on a GPU.

.. _python_spectral_install:

Spectral cones
""""""""""""""

To enable :ref:`spectral cone <spectral_cones>` support (log-determinant,
nuclear norm, :math:`\ell_1` norm, sum-of-largest-eigenvalues), install with:

.. code:: bash

  python -m pip install -Csetup-args=-Duse_spectral_cones=true .

This requires LAPACK (enabled by default). See
:ref:`Spectral cone keys <python_spectral_cone_keys>` for the cone dict
keys.

Testing
"""""""

To test that SCS installed correctly, and you have pytest installed, run

.. code:: bash

  python -m pytest .

See :ref:`here <python_interface>` for the full SCS python API.

OpenMP and the GPU indirect solver
----------------------------------

These are options of the Meson build, like the ones above. They are off by
default and the pre-built wheels are built without them, so they need an
install from source.

.. _python_openmp_install:

OpenMP and parallel PSD projection
""""""""""""""""""""""""""""""""""

To compile SCS with OpenMP, install with:

.. code:: bash

  python -m pip install -Csetup-args=-Duse_openmp=true .

SCS then uses several threads for its matrix-vector products and for the
projection onto exponential cones. The number of threads is the OpenMP thread
count: set the environment variable :code:`OMP_NUM_THREADS` before starting
Python.

To also project the positive semidefinite (PSD) blocks of a problem in
parallel, add the option :code:`use_parallel_psd`:

.. code:: bash

  python -m pip install -Csetup-args=-Duse_openmp=true -Csetup-args=-Duse_parallel_psd=true .

:code:`use_parallel_psd` requires :code:`use_openmp` and :code:`use_lapack`
(on by default). The build stops with an error without them. Read
:ref:`blas_lapack_openmp` before turning it on: the projection is only faster
when the BLAS runs on one thread inside the threads of SCS, and slower when it
does not. That section also says how to check that the blocks are projected in
parallel.

.. _python_gpu_indirect_install:

GPU indirect
""""""""""""

To build the :ref:`GPU indirect solver <gpu_indirect>`, with the CUDA toolkit
(including cuBLAS and cuSPARSE) installed, install with:

.. code:: bash

  python -m pip install -Csetup-args=-Duse_gpu=true -Csetup-args=-Dint32=true .

The :code:`int32=true` flag is required: the build stops with an error for a
GPU build with 64-bit integers. It applies to the whole build, so the CPU
solvers of that installation use 32-bit integers too. By default this build
does not store the transpose of :code:`A` in GPU memory. Add
:code:`-Csetup-args=-Dgpu_atrans=true` to store it, which uses more GPU memory
and makes the products with :code:`A` faster. Select the backend at runtime
with :code:`linear_solver=scs.LinearSolver.GPU_INDIRECT`.

The GPU indirect solver is effectively deprecated; the cuDSS direct solver
above is the recommended GPU backend.

