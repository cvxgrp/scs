.. _python_install:

Python
======

The easiest way to install the python version is using `pip <https://pypi.org/project/pip/>`_:

.. code:: bash

  pip install scs

SCS depends on python packages numpy and scipy to run and on numpy to install.
It uses numpy to tell it what :ref:`BLAS and LAPACK <blas_lapack>` libraries
to link against. If you run into an error like this:

.. code:: bash

  RuntimeError: Found /usr/lib/libcblas.dylib, but that file is a symbolic link to
  the MacOS Accelerate framework, which is not supported by NumPy

you can try:

.. code:: bash

  brew install openblas
  OPENBLAS="$(brew --prefix openblas)" pip install scs

You can also install directly from source

.. code:: bash

  git clone --recursive https://github.com/bodono/scs-python.git
  cd scs-python
  python setup.py install

You can install the MKL Pardiso interface using

.. code:: bash

  python setup.py install --scs --mkl

You can install the GPU interface using

.. code:: bash

  python setup.py install --scs --gpu

To test that SCS installed correctly, and you have pytest installed, run

.. code:: bash

  pytest

See :ref:`here <python_interface>` for the API.
