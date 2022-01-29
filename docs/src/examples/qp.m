% First, make sure SCS is in the path so MATLAB can call it
addpath("/Users/bodonoghue/git/scs-matlab")

% Set up data
data.P = sparse([3., -1.; -1., 2.]);
data.A = sparse([-1., 1.; 1., 0.; 0., 1.]);
data.b = [-1; 0.3; -0.5];
data.c = [-1.; -1.];
cone.z = 1;
cone.l = 2;

% Optional solver settings
settings = struct('eps_abs', 1e-9, 'eps_rel', 1e-9);

% Solve!
[x, y, s, info] = scs(data, cone, settings);

disp(sprintf("SCS took %i iters", info.iter))
disp("Optimal solution vector x*:");
disp(x)
disp("Optimal dual vector y*:");
disp(y)
