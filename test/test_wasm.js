import createSCS from '../out/scs.js';

const SCS = await createSCS();

const settings = new SCS.ScsSettings();
SCS.setDefaultSettings(settings);

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
        // bu: null,
        // bl: null,
        // bsize: 0,
        // q: null,
        // qsize: 0,
        // s: null,
        // ssize: 0,
        // ep: 0,
        // ed: 0,
        // p: null,
        // psize: 0
    };

    const settings = new SCS.ScsSettings();
    SCS.setDefaultSettings(settings);
    settings.epsAbs = 1e-9;
    settings.epsRel = 1e-9;

    // First solve without warm start
    const solution = SCS.solve(data, cone, settings);
    console.log('First solution:', solution);

    // Enable warm start
    settings.warmStart = true;
    settings.verbose = 1;

    // Second solve with warm start
    const solution2 = SCS.solve(data, cone, settings, solution);
    console.log('Second solution with warm start:', solution2);
});
