.. _matlab_install:

MATLAB
======

To install SCS in Matlab from source:

.. code:: bash

  git clone --recursive https://github.com/bodono/scs-matlab.git

Then in a Matlab session

.. code:: matlab

  cd <path/to/scs-matlab>
  make_scs

Remember to include the scs-matlab directory in your Matlab path if you wish to
use the mex file in your Matlab code.

.. code:: matlab

  addpath(pwd)
  savepath

See :ref:`here <matlab_interface>` for the API.

