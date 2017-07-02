rng('default');
rng(1);

n = 10000; m = ceil(n/5); s=ceil(n/10);

x_true=[randn(s,1);zeros(n-s,1)]; % true sparse signal
x_true = x_true(randperm(n));

density = 0.1;
rcA = .01;
A=sprandn(m,n,density,rcA);

b = A*x_true + 0.1*randn(m,1);
mu = 1;

%%
do_super_scs = 0;

tic;
cvx_begin
    cvx_solver scs
    cvx_solver_settings('eps', 1e-4,...
        'scale', 1,...
        'do_super_scs',do_super_scs,...
        'direction', 100,...
        'k0', 0,...
        'memory', 100,...
        'rho_x', 0.001,...
        'verbose', 2)
    variable x_c(n)
    minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
cvx_end
toc

%%
n=800;
P = randn(n,n);
tic;
cvx_begin sdp
cvx_solver scs
cvx_solver_settings('eps', 1e-4,...
        'scale', 1,...
        'do_super_scs',0,...
        'direction', 100,...
        'k0', 0,...
        'memory', 50,...
        'rho_x', 0.001,...
        'verbose', 2)
    variable Z(n,n) hermitian toeplitz
    dual variable Q
    minimize( norm( Z - P, 'fro' ) )
    Z >= 0 : Q;
cvx_end
toc
