.. _javascript_interface:

JavaScript / WebAssembly
========================

After :ref:`building the WebAssembly version <javascript_install>`, you can use SCS in JavaScript environments including browsers and Node.js.

Basic Usage
-----------

In Node.js, you can use SCS as a Node.js module:

.. code-block:: javascript

    const Module = require('scs.js');
    // or: import Module from 'scs.js';

    // Wait for WASM module initialization
    Module.onRuntimeInitialized = function() {
        // define problem here
        Module.solve(data, cone, settings);
    };

In browsers, you can load SCS using a script tag:

.. code-block:: html

    <script src="scs.js"></script>
    <script>
        Module.onRuntimeInitialized = function() {
            // define problem here
            Module.solve(data, cone, settings);
        };
    </script>

Data Format
-----------

Problem data must be provided as sparse matrices in CSC format using the following structure:

.. code-block:: javascript

    const data = {
        m: number,     // Number of rows of A
        n: number,     // Number of cols of A and of P
        A_x: number[], // Non-zero elements of matrix A
        A_i: number[], // Row indices of A elements
        A_p: number[], // Column pointers for A
        P_x: number[], // Non-zero elements of matrix P (optional)
        P_i: number[], // Row indices of P elements (optional)
        P_p: number[], // Column pointers for P (optional)
        b: number[],   // Length m array
        c: number[]    // Length n array
    };

One way to handle the CSC format in javascript is via the 
`Math.js library <https://mathjs.org/docs/reference/classes/sparsematrix.html>`_,
for example

.. code-block:: javascript

    // npm install mathjs
    const { matrix } = require('mathjs');
    // or import { matrix } from 'mathjs';
    // or <script src="https://unpkg.com/mathjs@14.0.1/lib/browser/math.js"></script>

    const A = matrix([
        [1, 0],
        [0, 1],
        [1, 1]
    ], 'sparse');

    const P = matrix([
        [3, 0],
        [0, 2]
    ], 'sparse');

    const data = {
        m: 3,
        n: 2,
        A_x: A._values,
        A_i: A._index,
        A_p: A._ptr,
        P_x: P._values,
        P_i: P._index,
        P_p: P._ptr,
        b: [-1.0, 0.3, -0.5],
        c: [-1.0, -1.0]
    };

Cone Specification
------------------

Cones are specified using the following structure:

.. code-block:: javascript

    const cone = {
        z: number,     // Number of linear equality constraints (primal zero, dual free)
        l: number,     // Number of positive orthant cones
        bu: number[],  // Upper box values (optional)
        bl: number[],  // Lower box values (optional)
        bsize: number, // Total length of box cone
        q: number[],   // Array of second-order cone constraints (optional)
        qsize: number, // Length of second-order cone array
        s: number[],   // Array of semidefinite cone constraints (optional)
        ssize: number, // Length of semidefinite constraints array
        ep: number,    // Number of primal exponential cone triples
        ed: number,    // Number of dual exponential cone triples
        p: number[],   // Array of power cone parameters (optional)
        psize: number  // Number of power cone triples convergence
    };

Settings
--------

Control solver behavior using settings:

.. code-block:: javascript

    const settings = new Module.ScsSettings();
    Module.setDefaultSettings(settings);

Available settings:

- ``normalize`` (boolean): Heuristically rescale problem data
- ``scale`` (number): Initial dual scaling factor
- ``adaptiveScale`` (boolean): Whether to adaptively update scale
- ``rhoX`` (number): Primal constraint scaling factor
- ``maxIters`` (number): Maximum iterations to take
- ``epsAbs`` (number): Absolute convergence tolerance
- ``epsRel`` (number): Relative convergence tolerance
- ``epsInfeas`` (number): Infeasible convergence tolerance
- ``alpha`` (number): Douglas-Rachford relaxation parameter
- ``timeLimitSecs`` (number): Time limit in seconds
- ``verbose`` (number): Output level (0-3)
- ``warmStart`` (boolean): Use warm starting

Solving Problems
----------------

Use the ``solve`` function to solve optimization problems:

.. code-block:: javascript

    const solution = Module.solve(data, cone, settings);

The solution object contains:

- ``x``: Primal variables
- ``y``: Dual variables
- ``s``: Slack variables
- ``info``: Solver information

    - ``iter``: Number of iterations
    - ``pobj``: Primal objective
    - ``dobj``: Dual objective
    - ``resPri``: Primal residual
    - ``resDual``: Dual residual
    - ``resInfeas``: Infeasibility residual
    - ``resUnbdd``: Unboundedness measure
    - ``solveTime``: Solve time
    - ``setupTime``: Setup time
- ``status``: Solution status code

Example
-------

Here's a complete example solving a quadratic program:

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