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

The installer automatically adds scs-matlab to your Matlab path and calls
``savepath`` so it persists across sessions. If ``savepath`` fails (e.g., due to
file permissions) it will print the ``addpath`` line you need to add to your
``startup.m``.

See :ref:`here <matlab_interface>` for the API.

