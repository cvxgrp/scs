rng('default');
rng(1);

n = 2000;
m = ceil(n/4);
density = 0.1;
G = sprandn(m,n,density);
f = randn(m,1) * n * density;
pow = 1.5;
tic
cvx_begin
    cvx_solver scs
    cvx_solver_settings('eps', 1e-4,...
        'verbose', 2,...
        'do_super_scs', 0, ...
        'rho_x', .001,  ...
        'k0', 0, ...
        'direction', 100, ...
        'memory', 50,...
        'use_indirect', 0);
    variable x(n)
    minimize(norm(x, pow))
    subject to
      G*x == f
    cvx_end
toc
