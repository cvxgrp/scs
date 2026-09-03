.. _python_install:

Python
======

The easiest way to install the python version is using `pip <https://pypi.org/project/pip/>`_:

.. code:: bash

  pip install scs

.. important::

   On x86-64 Linux the pre-built wheels include the :ref:`MKL Pardiso <mkl>`
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

- :code:`QDLDL` — the default sparse direct solver (bundled with SCS).
- :code:`CPU_INDIRECT` — the sparse matrix-free solver based on conjugate
  gradients.

The remaining backends require either a platform-specific library (Apple
Accelerate, MKL) or a build-time flag plus an external dependency (LAPACK for
dense, CUDA + cuDSS for GPU). Each section below describes what to install and
how to enable the backend. See :ref:`linear_solver` for an overview of each
solver and :ref:`python_interface` for how to select one at runtime.

Apple Accelerate (macOS)
""""""""""""""""""""""""

On macOS the Apple Accelerate backend is built and included automatically —
no extra install flags are needed. It uses the Accelerate framework's sparse
LDL\ :sup:`T` solver, which is optimized for Apple hardware including Apple
Silicon. The default :code:`linear_solver=scs.LinearSolver.AUTO` selects the
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

Build from source with:

.. code:: bash

  python -m pip install -Csetup-args=-Duse_lapack=true .

This requires BLAS and LAPACK development headers to be discoverable by
:code:`pkg-config`. Most platforms satisfy this out of the box (Apple
Accelerate on macOS, OpenBLAS / MKL on Linux, MKL on Windows). Select the
backend at runtime with
:code:`linear_solver=scs.LinearSolver.CPU_DENSE`.

GPU direct (cuDSS)
""""""""""""""""""

The :ref:`cuDSS backend <cudss_solver>` runs the sparse direct factorization
and solves on an NVIDIA GPU via `NVIDIA cuDSS
<https://developer.nvidia.com/cudss>`_. For large problems it is typically
substantially faster than any CPU backend.

**Prerequisites.** The build links against both the CUDA runtime and cuDSS,
so you need all of the following installed and discoverable by
:code:`pkg-config` / the linker before running :code:`pip install`:

1. An NVIDIA GPU with a recent CUDA-capable driver.
2. The `CUDA Toolkit <https://developer.nvidia.com/cuda-downloads>`_
   (provides :code:`nvcc`, the CUDA runtime headers, and :code:`cuda.pc`
   used by the build).
3. The `cuDSS library <https://developer.nvidia.com/cudss-downloads>`_ (ships
   a :code:`cudss.pc` pkg-config file). cuDSS is also available on
   `conda-forge <https://anaconda.org/conda-forge/libcudss>`_ as
   :code:`libcudss` / :code:`libcudss-dev`.

Make sure the directories containing :code:`cuda.pc` and :code:`cudss.pc` are
on :code:`PKG_CONFIG_PATH`, and that the corresponding shared libraries are on
:code:`LD_LIBRARY_PATH` (Linux) at runtime. A typical Linux environment looks
like:

.. code:: bash

  export PATH=/usr/local/cuda/bin:$PATH
  export PKG_CONFIG_PATH=/usr/local/cuda/lib64/pkgconfig:/opt/nvidia/cudss/lib64/pkgconfig:$PKG_CONFIG_PATH
  export LD_LIBRARY_PATH=/usr/local/cuda/lib64:/opt/nvidia/cudss/lib64:$LD_LIBRARY_PATH

**Install.** Once the prerequisites are in place, build SCS with:

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_cudss=true -Csetup-args=-Dint32=true .

The :code:`int32=true` flag is required because cuDSS only supports 32-bit
integer indices. Select the backend at runtime with
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
:ref:`python_spectral_cone_keys` for the cone dict keys.

Testing
"""""""

To test that SCS installed correctly, and you have pytest installed, run

.. code:: bash

  python -m pytest .

See :ref:`here <python_interface>` for the full SCS python API.

Legacy options
--------------

You can install with OpenMP parallelization support using

.. code:: bash

  python legacy_setup.py install --scs --openmp

You can install the :ref:`GPU indirect solver <gpu_indirect>` using

.. code:: bash

  python legacy_setup.py install --scs --gpu

The GPU indirect solver is effectively deprecated; the cuDSS direct solver
above is the recommended GPU backend.

