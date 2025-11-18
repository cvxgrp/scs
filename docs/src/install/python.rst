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

MKL
"""

If you have MKL, you can install the MKL Pardiso interface using

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_mkl=true .

See :ref:`here <python_interface>` for how to enable MKL when solving. MKL is
typically faster than the built-in linear system solver.

GPU
"""

If you have a GPU and cuDSS installed you can install the GPU direct sparse
solver using

.. code:: bash

  python -m pip install -Csetup-args=-Dlink_cudss=true -Csetup-args=-Dint32=true .

See :ref:`here <python_interface>` for how to enable the GPU when solving. The
sparse direct GPU solver is typically very fast.

See `here <https://colab.research.google.com/drive/1POCgDNFg8fycHMI9T9N6V3iHFhXRthjn?usp=sharing>`_ for an example colab where the cuDSS version of SCS, along with
required dependencies, is installed and used.

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

