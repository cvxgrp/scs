clc
clear all;

K = struct( 'f',0,...
            'l',0,...
            'q',4,...
            's',0,...
            'ep',0,...
            'ed',0,...
            'p',[] );

A(1,1) = 0.3;
A(4,1) = -0.5;
A(2,2) = 0.7;
A(4,2) = 0.9;
A(3,3) = 0.2;

b = [0.2; 0.1; -0.1; 0.1];
c = [1;-2;-3];

data.A = sparse(A);
data.b = b;
data.c = c;

params.eps=1e-3;
params.do_super_scs = 1;
params.max_iters = 2500;
params.alpha = 1.5;
params.scale = 1;
params.verbose = 1;
params.normalize=1;
params.direction = 5;
params.beta = 0.5;
params.c1 =1.0-1e-4;
params.c_bl = 0.999;
params.k0 = 1;
params.k1 = 1;
params.k2 = 1;
params.ls = 10;
params.sigma = 1e-2;
params.thetabar = 0.1;
params.rho_x = 1;
params.memory = 10;
params.sse = 0.7;
params.verbose = 2;

[x,y,s,info]=scs(data, K, params);