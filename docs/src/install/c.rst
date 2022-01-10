.. _c_install:

C / C++
=======

First clone the SCS repo from GitHub

.. code:: bash

  git clone https://github.com/cvxgrp/scs.git

CMake
^^^^^

Thanks to the `CMake <cmake.org>`_ buildsystem (contributed by `Giulio Romualdi
<https://github.com/GiulioRomualdi>`__), SCS can be easily compiled and linked
to other CMake projects. To use the cmake buld system please run the following
commands:

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
:code:`scs::scsindir` as well as a header file `scs.h` that defines the API. The
libraries can be imported using the find_package CMake command and used by
calling target_link_libraries as in the following example:

.. code:: bash

  cmake_minimum_required(VERSION 3.0)
  project(myproject)
  find_package(scs REQUIRED)
  add_executable(example example.cpp)

  # To use the direct method
  target_link_libraries(example scs::scsdir)

  # To use the indirect method
  target_link_libraries(example scs::scsindir)

Makefile
^^^^^^^^
Alternatively you can use the Makefile and manage the libraries yourself 

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

If you have a GPU and have CUDA installed, you can also execute make gpu to
compile SCS to run on the GPU which will create additional libraries and demo
binaries in the out folder corresponding to the GPU version.  Note that the GPU
(usually) requires 32 bit ints, which can be enforced by compiling with
:code:`DLONG=0`.

.. code:: bash

  make gpu DLONG=0
  out/run_tests_gpu_indirect

To use the libraries in your own source code, compile your code with the linker
option :code:`-L(PATH_TO_SCS_LIBS)` and :code:`-lscsdir` or :code:`-lscsindir`
(as needed). The API and required data structures are defined in the file
:code:`include/scs.h`. The four main API functions are:

