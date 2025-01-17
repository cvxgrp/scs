.. _javascript_example:

JavaScript / WebAssembly
========================

These examples assume that you have loaded `scs.js`, either in Node.js via

.. code-block:: javascript

    const createSCS = require('scs-solver'); // if using CommonJS
    import createSCS from 'scs-solver'; // if using ES6 modules

or in the browser via a script tag or an ES6 module import (see :ref:`the 
install page <javascript_install>`).

Live Demo
---------

Here is a live demo, which computes the point :math:`x` that is closest to the 
origin, subject to lying in the half-plane :math:`x_1 + x_2 \ge b` and within 
a disc of radius :math:`r` centered at :math:`c = (-2, 2.5)`. Here, the values 
of :math:`b` and :math:`r` are set by sliders. The solution is computed using 
two second-order cones, one for the disc constraint, and one to encode
:math:`\|x\|_2 \leq d`, where :math:`d` is the distance from the origin, which
is the quantity to be minimized.

.. raw:: html
   :file: javascript_disc.html

.. raw:: html

  <details style="margin: 1em 0;">
    <summary style="cursor: pointer; margin-bottom: 1em;">Show code</summary>

.. literalinclude:: javascript_disc.html
  :language: html

.. raw:: html

  </details>

Basic Usage
-----------

Here's the :ref:`basic example from C <c_example>` translated to JavaScript:

.. code-block:: javascript

    createSCS().then(SCS => {
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

        const settings = new SCS.ScsSettings();
        SCS.setDefaultSettings(settings);
        settings.epsAbs = 1e-9;
        settings.epsRel = 1e-9;

        const solution = SCS.solve(data, cone, settings);
        console.log(solution);

        // re-solve using warm start (will be faster)
        settings.warmStart = true;
        const solution2 = SCS.solve(data, cone, settings, solution);
    });

This prints the solution object to the console:

.. code-block:: javascript

    {
      x: [ 0.3000000000043908, -0.6999999999956144 ],
      y: [ 2.699999999995767, 2.0999999999869825, 0 ],
      s: [ 0, 0, 0.1999999999956145 ],
      info: {
        iter: 100,
        status: 'solved',
        linSysSolver: 'sparse-direct-amd-qdldl',
        statusVal: 1,
        scaleUpdates: 0,
        pobj: 1.2349999999907928,
        dobj: 1.2350000000001042,
        resPri: 4.390808429506794e-12,
        resDual: 1.4869081633461182e-13,
        gap: 9.311465734712679e-12,
        resInfeas: 1.3043478260851176,
        resUnbddA: NaN,
        resUnbddP: NaN,
        compSlack: 0,
        setupTime: 2.796667,
        solveTime: 0.505584,
        scale: 0.1,
        rejectedAccelSteps: 0,
        acceptedAccelSteps: 0,
        linSysTime: 0.047704000000000024,
        coneTime: 0.07804600000000002,
        accelTime: 0
      },
      statusVal: 1,
      status: 'SOLVED'
    }

Entropy Example
---------------

Next, we will consider a problem involving maximum entropy. Given a vector 
:math:`y \in \mathbf{R}^n`, we want to optimize a function involving entropy
over the unit simplex.

.. math::
  \begin{align*}
	  \text{minimize} \quad & \sum_{i = 1}^n x_i \log x_i - \langle y, x \rangle \\
    \text{subject to} \quad & \sum_{i = 1}^n x_i = 1 \\
    & x \geq 0
  \end{align*}

It is known that for the optimal solution, we have :math:`x_i \propto e^{y_i}`.

This problem can be formulated using the :ref:`(primal) exponential cone <cones>`,
defined as 

.. math::
  \begin{align*}
    \mathcal{K}_{\text{exp}} &= \{ (x,y,z) \in \mathbf{R}^3 \mid y e^{x/y} \leq z, y>0  \} \\
    &= \{ (x,y,z) \in \mathbf{R}^3 \mid y \log(z/y) \geq x, y>0, z>0 \}
  \end{align*}

Our formulation is then:

.. math::
  \begin{align*}
    \text{minimize} \quad & \sum_{i = 1}^n t_i - \langle y, x \rangle \\
    \text{subject to} \quad & \sum_{i = 1}^n x_i = 1 \\
    & x_i \geq 0 \: && \forall i \\
    & (-t_i, x_i, 1) \in \mathcal{K}_{\text{exp}} \: && \forall i
  \end{align*}

To implement this problem in JavaScript, we will use the sparse matrix
implementation from the `Math.js library <https://mathjs.org/docs/reference/classes/sparsematrix.html>`_.

.. code-block:: javascript

    const createSCS = require('./out/scs.js');
    const math = require('./math.js');

    createSCS().then(SCS => {
        const n = 5;
        const y = Array.from({ length: n }, () => Math.random());

        const A = math.matrix('sparse');
        const b = [];

        let constraintIndex = 0;

        const x_vars = Array.from({ length: n }, (_, i) => i);
        const t_vars = Array.from({ length: n }, (_, i) => i + n);

        // equality constraint (zero cone)
        let numEqCones = 0;
        for (let i = 0; i < n; i++) {
            A.set([constraintIndex, x_vars[i]], 1);
        }
        b.push(1);
        constraintIndex++;
        numEqCones++;

        // inequality constraints (positive cone)
        let numPosCones = 0;
        for (let i = 0; i < n; i++) {
            A.set([constraintIndex, x_vars[i]], -1);
            b.push(0);
            constraintIndex++;
            numPosCones++;
        }

        // exponential cone constraints
        let numExpCones = 0;
        for (let i = 0; i < n; i++) {
            // (-t_i, x_i, 1) in exponential cone
            A.set([constraintIndex, t_vars[i]], 1);
            b.push(0);
            constraintIndex++;
            A.set([constraintIndex, x_vars[i]], -1);
            b.push(0);
            constraintIndex++;
            // last element is constant, so A has a 0-row; set arbitrary index to 0
            A.set([constraintIndex, x_vars[i]], 0);
            b.push(1);
            constraintIndex++;
            numExpCones++;
        }

        // objective function
        const c = Array.from({ length: 2 * n }, (_, i) => 0);
        for (let i = 0; i < n; i++) {
            c[x_vars[i]] = -y[i];
            c[t_vars[i]] = 1;
        }

        const data = {
            m: A._size[0],
            n: A._size[1],
            A_x: A._values,
            A_i: A._index,
            A_p: A._ptr,
            b: b,
            c: c,
        };

        const cone = {
            z: numEqCones,
            l: numPosCones,
            ep: numExpCones,
        };

        const settings = new SCS.ScsSettings();
        SCS.setDefaultSettings(settings);
        settings.epsAbs = 1e-9;
        settings.epsRel = 1e-9;

        const solution = SCS.solve(data, cone, settings);
        console.log("SCS solution:", solution.x.slice(0, n));

        const denominator = y.map(y_i => Math.exp(y_i)).reduce((a, b) => a + b, 0);
        const predicted_solution = y.map(y_i => Math.exp(y_i) / denominator);
        console.log("Predicted solution:", predicted_solution);
    });