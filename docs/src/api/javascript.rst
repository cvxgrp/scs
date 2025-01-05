.. _javascript_interface:

JavaScript / WebAssembly
========================

After :ref:`building the WebAssembly version <javascript_install>`, you can use SCS in JavaScript environments including browsers and Node.js.

Basic Usage
----------

In Node.js, you can use SCS as a Node.js module:

.. code-block:: javascript

    const Module = require('scs.js');

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
----------

Problem data must be provided as sparse matrices in CSC format using the following structure:

.. code-block:: javascript

    const data = {
        m: number,     // Number of constraints
        n: number,     // Number of variables
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
[Math.js library](https://mathjs.org/docs/reference/classes/sparsematrix.html),
for example

.. code-block:: javascript

    const math = require('mathjs');

    const A = math.sparse([
        [1, 0],
        [0, 1],
        [1, 1]
    ]);

    const P = math.sparse([
        [3, 0],
        [0, 2]
    ]);

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
----------------

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
--------------

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
            bsize: 0,
            qsize: 0,
            ssize: 0,
            ep: 0,
            ed: 0,
            psize: 0
        };

        const settings = new Module.ScsSettings();
        Module.setDefaultSettings(settings);
        settings.epsAbs = 1e-9;
        settings.epsRel = 1e-9;

        const solution = Module.solve(data, cone, settings);
        console.log(solution);
    };