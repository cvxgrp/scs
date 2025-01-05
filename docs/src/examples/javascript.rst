.. _javascript_example:

JavaScript / WebAssembly
========================

These examples assume that you have loaded `scs.js`, either in Node.js via

.. code-block:: javascript

    const Module = require('scs.js'); // if using CommonJS
    import Module from 'scs.js'; // if using ES6 modules

of in the browser via

.. code-block:: html

    <script src="scs.js"></script>

Basic Usage
-----------

Here's the :ref:`basic example from C <c_example>` translated to JavaScript:

.. code-block:: javascript

    Module.onRuntimeInitialized = function() {
        const data = {
            m: 3,
            n: 2,
            A_x: [-1.0, 1.0, 1.0, 1.0],
            A_i: [0, 1, 0, 2],
            A_p: [0, 2, 4],
            P_x: [3.0, -1.0, 2.0],
            P_i: [0, 0, 1],
            P_p: [0, 1, 3],
            b: [-1.0, 0.3, -0.5],
            c: [-1.0, -1.0]
        };

        const cone = {
            z: 1,
            l: 2,
        };

        const settings = new Module.ScsSettings();
        Module.setDefaultSettings(settings);
        settings.epsAbs = 1e-9;
        settings.epsRel = 1e-9;

        const solution = Module.solve(data, cone, settings);
        console.log(solution);

        // re-solve using warm start (will be faster)
        settings.warmStart = true;
        const solution2 = Module.solve(data, cone, settings, solution);
    };

This prints the solution object to the console:

.. code-block:: javascript

    {
      x: [ 0.3000000000043908, -0.6999999999956144 ],
      y: [ 2.699999999995767, 2.0999999999869825, 0 ],
      s: [ 0, 0, 0.1999999999956145 ],
      info: {
        iter: 100,
        pobj: 1.2349999999907928,
        dobj: 1.2350000000001042,
        resPri: 4.390808429506794e-12,
        resDual: 1.4869081633461182e-13,
        resInfeas: 1.3043478260851176,
        resUnbdd: NaN,
        solveTime: 0.598459,
        setupTime: 11.603125
      },
      status: 1
    }