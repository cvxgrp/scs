rng('default');
rng(1);

m = 3000; 
n = ceil(m/2);
A = sprandn(m,n,0.5,0.01);
b = 10*randn(m,1);
G = 2*sprandn(2*n, n, 0.1);



%%
tic;
cvx_begin
    cvx_solver scs
    cvx_solver_settings('eps', 1e-3,...
        'verbose', 1,...
        'do_super_scs', 1, ...
        'rho_x', .001,  ...
        'k0', 0, ...
        'direction', 100, ...
        'memory', 20);
    variable x(n)
    minimize( norm(A*x-b) )
    subject to
        norm(G*x) <= 1
cvx_end
toc
