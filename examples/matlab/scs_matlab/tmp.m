clear all;close all;
randn('seed',0)

n= 100;

A = randn(n); A = A*A*A;
b = randn(n,1);
c = randn(n,1);
lam = 1;

%%

cvx_begin
cvx_solver scs_matlab
%cvx_solver('scs')
cvx_solver_settings('anderson_lookback', 30, 'extra_verbose', 0, 'eps',1e-9, 'scale', 1, 'normalize', 1, 'use_indirect', 0, 'gen_plots', 1, 'max_iters', 100000)
variable x(n)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end

%%

cvx_begin
cvx_solver scs_acc
%cvx_solver('scs')
cvx_solver_settings('anderson_lookback', 10, 'extra_verbose', 0, 'eps',1e-9, 'scale', 1, 'normalize', 1, 'use_indirect', 1, 'gen_plots', 1, 'max_iters', 100000)
variable x(n)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end

%%
cvx_begin
variable x(n)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end


