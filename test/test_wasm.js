// Import the WASM module
const Module = require('../out/scs.js');

// Wait for the module to be initialized
Module.onRuntimeInitialized = async function() {
    // Define problem data (same as demo.html)
    const P_x = [3.0, -1.0, 2.0];
    const P_i = [0, 0, 1];
    const P_p = [0, 1, 3];
    const A_x = [-1.0, 1.0, 1.0, 1.0];
    const A_i = [0, 1, 0, 2];
    const A_p = [0, 2, 4];
    const b = [-1.0, 0.3, -0.5];
    const c = [-1.0, -1.0];

    // Create ScsData object
    const data = {
        m: 3,
        n: 2,
        A_x, A_i, A_p,
        P_x, P_i, P_p,
        b, c
    };

    // Create ScsCone object
    const cone = {
        z: 1,
        l: 2,
        bu: null,
        bl: null,
        bsize: 0,
        q: null,
        qsize: 0,
        s: null,
        ssize: 0,
        ep: 0,
        ed: 0,
        p: null,
        psize: 0
    };

    // Initialize settings
    const settings = new Module.ScsSettings();
    Module.setDefaultSettings(settings);
    settings.epsAbs = 1e-9;
    settings.epsRel = 1e-9;
    settings.verbose = 1;  // Set to 1 to see solver output

    // Solve the problem
    const solution = Module.solve(data, cone, settings);

    // Print results
    console.log('Solution status:', solution.status);
    console.log('\nPrimal variables (x):');
    console.log(Array.from(solution.x));
    console.log('\nDual variables (y):');
    console.log(Array.from(solution.y));
    console.log('\nSlack variables (s):');
    console.log(Array.from(solution.s));
    console.log('\nSolver info:');
    console.log(solution.info);
};
