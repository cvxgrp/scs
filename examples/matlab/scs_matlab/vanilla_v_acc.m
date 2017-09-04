clear all; close all
K = struct('f',50,'l',300,'q',[10],'s',[20],'ep',0,'ed',0,'p',[]);

density = 0.1; % A matrix density

m = getConeDims(K);
n = round(m/3);
params = struct('use_indirect', 1, 'eps', 1e-9, 'normalize', 1, 'scale', 1, 'cg_rate',2, 'max_iters', 2500, 'alpha', 1.5, 'gen_plots', 1);

%% generate primal-dual feasible cone prob:
% Ax + s = b, s \in K, A'y + c = 0, y \in K*, s'y = 0
z = randn(m,1);
y = proj_dual_cone(z,K); % y = s - z;
s = y - z; %s = proj_cone(z,K);

A = sprandn(m,n,density);
[U,S,V] = svds(A);
mS = min(diag(S));
S = mS*((S/mS).^70);
A = sparse(U*S*V');
x = randn(n,1);
c = -A'*y;
b = A*x + s;

data.A = A;
data.b = b;
data.c = c;

params.anderson_lookback=10;
[xi,yi,si,infoi] = scs_acc(data,K,params);
c'*x
(c'*xi - c'*x) / (c'*x)
b'*y
(b'*yi - b'*y) / (b'*y)

rmpath ~/Dropbox/research/superscs/matlab
addpath ~/Dropbox/research/scs/matlab
[xi,yi,si,infoi] = scs(data,K,params);
c'*x
(c'*xi - c'*x) / (c'*x)
b'*y
(b'*yi - b'*y) / (b'*y)

addpath ~/Dropbox/research/superscs/matlab
rmpath ~/Dropbox/research/scs/matlab
[xi,yi,si,infoi] = scs(data,K,params);
c'*x
(c'*xi - c'*x) / (c'*x)
b'*y
(b'*yi - b'*y) / (b'*y)
