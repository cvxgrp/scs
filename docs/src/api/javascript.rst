.. _javascript_interface:

JavaScript / WebAssembly
========================

After :ref:`installing the scs-solver package <javascript_install>`, or after 
loading ``scs.js`` :ref:`in the browser <javascript_install>`, you can use SCS
using the following interface.

Instantiating the Solver
------------------------

In Node.js, you can use SCS as follows using CommonJS:

.. code-block:: javascript

    const createSCS = require('scs-solver');

    createSCS().then(SCS => {
        // define problem here
        SCS.solve(data, cone, settings);
    });

Alternatively, you can use ES6 modules, as well as async/await:

.. code-block:: javascript

    import createSCS from 'scs-solver';

    async function main() {
        const SCS = await createSCS();
        // define problem here
        SCS.solve(data, cone, settings);
    }

    main();

In browsers, you can load SCS using a script tag:

.. code-block:: html

    <script src="https://unpkg.com/scs-solver/dist/scs.js"></script>
    <script>
        createSCS().then(SCS => {
            // define problem here
            SCS.solve(data, cone, settings);
        });
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
        z: number,     // Number of zero cones
        l: number,     // Number of positive (or linear) cones
        bu: number[],  // Box cone upper values
        bl: number[],  // Box cone lower values
        bsize: number, // Total length of box cone
        q: number[],   // Array of second-order cone lengths
        qsize: number, // Number of second-order cones
        ep: number,    // Number of primal exponential cone triples
        ed: number,    // Number of dual exponential cone triples
        p: number[],   // Array of power cone parameters
        psize: number  // Number of power cone triples
    };

Note that positive semidefinite cones are not supported in the JavaScript interface.

Usually, not all cone types are used in a problem, in which case the unused 
cones can be omitted. For example, if only zero and positive cones are used:

.. code-block:: javascript

    const cone = {
        z: 1,
        l: 2
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

    const solution = Module.solve(data, cone, settings, [warmStartSolution]);

The function takes an optional ``warmStartSolution`` object to warm-start the solver,
provided ``settings.warmStart`` is set to ``true``.

The returned ``solution`` object contains:

- ``x``: Primal variables
- ``y``: Dual variables
- ``s``: Slack variables
- ``info``: Solver information

    - ``iter``: Number of iterations
    - ``pobj``: Primal objective
    - ``dobj``: Dual objective
    - ``solveTime``: Solve time
    - and :ref:`other solver information <info>`
    
- ``status``: Solution status (e.g. ``SOLVED``, ``INFEASIBLE``, ``UNBOUNDED``, see :ref:`exit flags <exit_flags>`)
- ``statusVal``: Solution status value (see :ref:`exit flags <exit_flags>`)

