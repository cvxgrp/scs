clear all;close all;
randn('seed', 81887)
rand('seed', 77)

n = 1000;

A = randn(n);
b0 = sprandn(n, 1, 0.3);
b = A * b0 + 0.05 * randn(n,1);
%b = randn(n,1);
lam = 1;

use_indirect = 1;

%%
%{%
cvx_begin
cvx_solver scs_matlab
cvx_solver_settings('extra_verbose', 0, 'eps',1e-9, 'scale', 1, 'normalize', 1, 'use_indirect', use_indirect, 'gen_plots', 0, 'max_iters', 5000)
variable x(n)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end
%}
%%
%{%
cvx_begin
cvx_solver scs_acc
%cvx_solver('scs_matlab')
cvx_solver_settings('anderson_lookback', 20, 'extra_verbose', 0, 'eps',1e-9, 'scale', 1, 'normalize', 1, 'use_indirect', use_indirect, 'gen_plots', 0, 'max_iters', 5000)
variable x(n)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end
%}
%%
rmpath ~/Dropbox/research/superscs/matlab
rmpath ~/Dropbox/research/scs/matlab
addpath ~/Dropbox/research/scs_clean/matlab
cvx_begin
variable x(n)
cvx_solver('scs')
cvx_solver_settings('extra_verbose', 0, 'eps',1e-7, 'scale', 1, 'normalize', 1, 'use_indirect', use_indirect, 'gen_plots', 1, 'max_iters', 10000)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end


%%
rmpath ~/Dropbox/research/superscs/matlab
rmpath ~/Dropbox/research/scs_clean/matlab
addpath ~/Dropbox/research/scs/matlab
cvx_begin
variable x(n)
cvx_solver('scs')
cvx_solver_settings('extra_verbose', 0, 'eps',1e-7, 'scale', 1, 'normalize', 1, 'use_indirect', use_indirect, 'gen_plots', 0, 'max_iters', 10000)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end

%%

rmpath ~/Dropbox/research/scs/matlab
rmpath ~/Dropbox/research/scs_clean/matlab
addpath ~/Dropbox/research/superscs/matlab
cvx_begin
variable x(n)
cvx_solver('scs')
cvx_solver_settings('extra_verbose', 0, 'eps',1e-7, 'scale', 1, 'normalize', 1, 'use_indirect', use_indirect, 'gen_plots', 1, 'max_iters', 10000)
minimize(norm(A*x - b) + lam*norm(x,1))
x >= 0
cvx_end



