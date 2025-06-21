.. _c_install:

C / C++
=======

First clone the SCS repo from GitHub

.. code:: bash

  git clone https://github.com/cvxgrp/scs.git

CMake
^^^^^

Thanks to the `CMake <https://cmake.org/cmake/help/latest/>`__ buildsystem SCS
can be easily compiled and linked to other CMake projects. To use the cmake
build system please run the following commands:

.. code:: bash

  cd scs
  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX:PATH=<custom-folder> ../
  make
  make install

You may also want to compile the tests. In this case when you configure the
project, please call the following command

.. code:: bash

  cmake -DCMAKE_INSTALL_PREFIX:PATH=<custom-folder> -DBUILD_TESTING=ON ../
  make
  ctest

Some :ref:`compile flags <compile_flags>` can be overridden using the
command line, for example we can compile the library (and headers) to use 64 bit
integers using:

.. code:: bash

  cmake -DCMAKE_INSTALL_PREFIX:PATH=<custom-folder> -DDLONG=ON ../
  make

By default the build-system will compile the library as shared. If you want to
compile it as static, please call the following command when you configure the
project

.. code:: bash

  cmake -DCMAKE_INSTALL_PREFIX:PATH=<custom-folder> -BUILD_SHARED_LIBS=OFF ../
  make

The CMake build-system exports two CMake targets called :code:`scs::scsdir` and
:code:`scs::scsindir` as well as a header file :code:`scs.h` that defines the
API.

If `MKL
<https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html>`_
is installed in your system and the :code:`MKLROOT` environment variable is
set, then additionally CMake will build and install the :ref:`MKL Pardiso
<mkl>` linear solver with target :code:`scs::scsmkl`.  (Note that the choice of
MKL compiler flags might not be right for your system and may need to be
`modified
<https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html>`_).


If you have a GPU and CUDA toolkit installed, along with the
`cuDSS <https://developer.nvidia.com/cudss>`_ library, you can compile SCS
with cuDSS support using CMake. First, ensure that the :code:`CUDA_PATH` and
:code:`CUDSS_PATH` environment variables are set, then configure with cuDSS enabled:

.. code:: bash

  cmake -DCMAKE_INSTALL_PREFIX:PATH=<custom-folder> -DUSE_CUDSS=ON -DDLONG=OFF ../
  make

Currently cuDSS only supports 32 bit integers (for sparse matrix idicies) so
:code:`DDLONG=OFF` is mandatory.
This will build and install the cuDSS linear solver with target :code:`scs::scscudss`.

The libraries can be imported using the find_package CMake command and used
by calling target_link_libraries as in the following example:

.. code:: bash

  cmake_minimum_required(VERSION 3.0)
  project(myproject)
  find_package(scs REQUIRED)
  add_executable(example example.cpp)

  # To use the direct method
  target_link_libraries(example scs::scsdir)

  # To use the indirect method
  target_link_libraries(example scs::scsindir)

  # To use the MKL Pardiso direct method
  target_link_libraries(example scs::scsmkl)

  # To use the cuDSS direct method
  target_link_libraries(example scs::scscudss)

Makefile
^^^^^^^^
Alternatively you can use the Makefile and manage the libraries and header
files yourself.  The public header files are :code:`scs.h` and
:code:`scs_types.h`.

.. code:: bash

  cd scs
  make


To compile and run the tests execute

.. code:: bash

  make test
  out/run_tests_direct
  out/run_tests_indirect

If make completes successfully, it will produce two static library files,
:code:`libscsdir.a`, :code:`libscsindir.a`, and two dynamic library files
:code:`libscsdir.ext`, :code:`libscsindir.ext` (where :code:`.ext` extension is
platform dependent) in the :code:`out` folder.

If `MKL
<https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html>`_
is installed in your system and the :code:`MKLROOT` environment variable is set,
then you can compile and test the :ref:`MKL Pardiso <mkl>` version of SCS using:

.. code:: bash

  make mkl
  out/run_tests_mkl

This will produce static library :code:`libscsmkl.a` and dynamic library
:code:`libscsmkl.ext` (again :code:`.ext` is platform dependent) in the
:code:`out` folder.  (Note that the choice of MKL compiler flags might not be right
for your system and may need to be `modified
<https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html>`_).

If you have a GPU and have CUDA installed, you can also execute make gpu to
compile SCS to run on the GPU which will create additional libraries and demo
binaries in the out folder corresponding to the GPU version.  Note that the GPU
(usually) requires 32 bit ints, which can be enforced by compiling with
:code:`DLONG=0`.

.. code:: bash

  make gpu DLONG=0
  out/run_tests_gpu_indirect

Finally, to compile and test the :ref:`cuDSS solver <cudss>` you need to have
CUDA toolkit, the :code:`nvcc` compiler, and
`cuDSS <https://developer.nvidia.com/cudss>`_ library installed.
Then set :code:`CUDA_PATH` and :code:`CUDSS_PATH` and execute

.. code:: bash

  make cudss DLONG=0
  out/run_tests_cudss

Currently cuDSS only supports 32 bit integers (for sparse matrix idicies) so
:code:`DLONG=0` is mandatory (see `the docs of cuDSS CSR matrix
<https://docs.nvidia.com/cuda/cudss/functions.html#cudssmatrixcreatecsr>`_).

To use the libraries in your own source code, compile your code with the linker
option :code:`-L(PATH_TO_SCS_LIBS)` and :code:`-lscsdir` or :code:`-lscsindir`
(as needed). The API and required data structures are defined in the file
:code:`include/scs.h` and documented :ref:`here <c_interface>`.

