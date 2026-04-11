.. _matlab_install:

MATLAB
======

The easiest way to install SCS in MATLAB is via the pre-compiled toolbox file.

Download Toolbox (.mltbx) - Recommended
---------------------------------------

The pre-compiled MATLAB Toolbox file includes binaries for Windows, Linux, and Apple Silicon Macs.

1. Go to the `SCS MATLAB interface releases <https://github.com/bodono/scs-matlab/releases>`_ page.
2. Download the latest ``SCS.mltbx`` file.
3. Open the file in MATLAB (or double-click it) to install.

Build from Source
-----------------

If you are on an unsupported platform or prefer to build from source:

1. Clone the repository recursively:

.. code:: bash

  git clone --recursive https://github.com/bodono/scs-matlab.git

2. In a MATLAB session:

.. code:: matlab

  cd <path/to/scs-matlab>
  make_scs

The installer automatically adds scs-matlab to your MATLAB path and calls
``savepath`` so it persists across sessions. If ``savepath`` fails (e.g., due to
file permissions) it will print the ``addpath`` line you need to add to your
``startup.m``.

See :ref:`here <matlab_interface>` for the API.
