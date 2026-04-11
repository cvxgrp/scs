.. _python_install:

Python
======

The easiest way to install the python version is using `pip <https://pypi.org/project/pip/>`_:

.. code:: bash

  pip install scs

You can also install directly from source

.. code:: bash

  git clone --recursive https://github.com/bodono/scs-python.git
  cd scs-python
  python -m pip install .

Apple Accelerate (macOS)
""""""""""""""""""""""""

On macOS the Apple Accelerate backend is built and included automatically —
no extra install flags are needed. It uses the Accelerate framework's sparse
LDL\ :sup:`T` solver, which is optimized for Apple hardware including Apple
Silicon. When using the default :code:`linear_solver=scs.LinearSolver.AUTO`,
Accelerate is selected automatically on macOS.

MKL
"""

If you have MKL, you can install the MKL Pardiso interface using

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_mkl=true .

When using the default :code:`linear_solver=scs.LinearSolver.AUTO`, MKL is
selected automatically on Linux and Windows if installed. MKL is
typically faster than the built-in QDLDL linear system solver.

GPU
"""

If you have a GPU and cuDSS installed you can install the GPU direct sparse
solver using

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_cudss=true -Csetup-args=-Dint32=true .

Select it at runtime with :code:`linear_solver=scs.LinearSolver.CUDSS`. The
sparse direct GPU solver is typically very fast.

See `here <https://colab.research.google.com/drive/1POCgDNFg8fycHMI9T9N6V3iHFhXRthjn?usp=sharing>`_ for an example colab where the cuDSS version of SCS, along with
required dependencies, is installed and used.

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

You can install the GPU indirect solver using

.. code:: bash

  python legacy_setup.py install --scs --gpu

