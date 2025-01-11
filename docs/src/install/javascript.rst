.. _javascript_install:

JavaScript / WebAssembly
========================

SCS is available in `a WebAssembly version 
<https://github.com/DominikPeters/scs.wasm>`_ that can be used with JavaScript.
Note that the JavaScript version does not support compiling with BLAS 
and LAPACK, so it does not support solving SDPs.

For building the WebAssembly version from source, see the `scs.wasm 
<https://github.com/DominikPeters/scs.wasm>`_ repository.

The page on the :ref:`JS interface <javascript_interface>` provides more 
information on how to use SCS in JavaScript environments once loaded or installed.

Install with npm
----------------

The package can be installed using `npm <https://www.npmjs.com/package/scs-solver>`_:

.. code:: bash

  npm install scs-solver

It can be used in Node.js and in the browser.

Loading from a CDN in the browser
---------------------------------

The package can also directly be included in a webpage using a CDN, by
using one of the following script tags:

.. code:: html

  <script src="https://unpkg.com/scs-solver/dist/scs.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/scs-solver/dist/scs.js"></script>

It can also be imported in JavaScript code using ES6 modules:

.. code:: html

  <script type="module">
    import createSCS from 'https://unpkg.com/scs-solver/dist/scs.mjs';
    // ...
  </script>

You can also host the files yourself, by `downloading the files 
<https://unpkg.com/browse/scs-solver/dist/>`_, and putting ``scs.js`` (or 
``scs.mjs``) and ``scs.wasm`` in the same directory.