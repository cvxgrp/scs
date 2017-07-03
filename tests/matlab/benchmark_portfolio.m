clear; clc;

% Generate portfolio optimization problem
rng('default');
rng(1);

density = 0.1;
rc = 0.5; % estimated reciprocal condition number

n = 100000;
m = 100;

mu = exp(0.01 * randn(n, 1)) - 1; % returns
D = rand(n,1) / 10; % idiosyncratic risk
F = sprandn(n, m, density, rc) / 10; % factor model
gamma = 1;
B = 1;

%% 
do_super_scs = 1;

tic;
cvx_begin
    cvx_solver scs
    cvx_solver_settings('eps', 1e-3,...
        'scale', 1,...
        'do_super_scs',do_super_scs,...
        'direction', 100,...
        'k0', 0,...
        'ls', 5,...
        'memory', 200,...
        'rho_x', 0.1,...
        'verbose', 2)
    variable x(n)
    maximize(mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
    sum(x) == B
    x >= 0
cvx_end
toc;
