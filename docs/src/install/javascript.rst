.. _javascript_install:

JavaScript / WebAssembly
========================

To build the WebAssembly (WASM) version of SCS, you will need to have `Emscripten` installed, 
and `emcc` in your path. To install it, run

.. code:: bash

  git clone https://github.com/emscripten-core/emsdk.git
  cd emsdk
  ./emsdk install latest
  ./emsdk activate latest
  source ./emsdk_env.sh

Now clone the SCS repo from GitHub and compile the WebAssembly (WASM) version of SCS.

.. code:: bash

  git clone https://github.com/cvxgrp/scs.git
  cd scs
  make wasm

If make completes successfully, you will find the compiled `scs.wasm` file 
and the JavaScript wrapper `scs.js` in the `out` directory. You can now use
these files in your JavaScript project, either in the browser or in Node.js.
