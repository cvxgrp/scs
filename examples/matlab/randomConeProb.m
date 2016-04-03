clear all;close all;

addpath('../../matlab')
addpath('./scs_matlab')

%cd '../../matlab'; make_scs; cd '../examples/matlab';

randn('seed',0);rand('seed',0);

%% generate random cone problem (solution NOT necessariy unique):

%%% types of probs to solve:
gen_feasible = true;
gen_infeasible = true;
gen_unbounded = true;
%%% solvers to test:
run_indirect = true;
run_direct = true;
run_cvx = false; % won't work if ep or ed > 0
cvx_solver = 'sdpt3';
run_scs_matlab = false; % SCS implemented in MATLAB

% set cone sizes (ep = ed = 0 if you want to compare against cvx):
K = struct('f',50,'l',60,'q',[0;1;0;2;100;20],'s',[0;2;1;2;20],'ep',10,'ed',15,'p',[-0.3, 0.25, 0.75, -0.4]);

density = 0.1; % A matrix density

m = getConeDims(K);
n = round(m/3);
params = struct('eps', 1e-3, 'normalize', 1, 'scale', 5, 'cg_rate',2, 'max_iters', 2500, 'alpha', 1.5);

%% generate primal-dual feasible cone prob:
% Ax + s = b, s \in K, A'y + c = 0, y \in K*, s'y = 0
if (gen_feasible)
    z = randn(m,1);
    y = proj_dual_cone(z,K); % y = s - z;
    s = y - z; %s = proj_cone(z,K);
    
    A = sprandn(m,n,density);
    x = randn(n,1);
    c = -A'*y;
    b = A*x + s;
    
    data.A = A;
    data.b = b;
    data.c = c;
    
    %cd '../../matlab'; write_scs_data(data,K,params,'randomConeFeasible'); cd '../examples/matlab';
    
    %indirect
    if (run_indirect)
        [xi,yi,si,infoi] = scs_indirect(data,K,params);
        c'*x
        (c'*xi - c'*x) / (c'*x)
        b'*y
        (b'*yi - b'*y) / (b'*y)
    end
    if (run_direct)
        % direct:
        [xd,yd,sd,infod] = scs_direct(data,K,params);
        c'*x
        (c'*xd - c'*x) / (c'*x)
        b'*y
        (b'*yd - b'*y) / (b'*y)
    end
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
    if (run_scs_matlab)
        params.use_indirect = true;
        [xi_m,yi_m,si_m,infoi_m] = scs_matlab(data,K,params);
        c'*x
        (c'*xi_m - c'*x) / (c'*x)
        b'*yi_m
        (b'*yi_m - b'*y) / (b'*y)
        params.use_indirect = false;
        [xd_m,yd_m,sd_m,infod_m] = scs_matlab(data,K,params);
        c'*x
        (c'*xd_m - c'*x) / (c'*x)
        b'*y
        (b'*yd_m - b'*y) / (b'*y)
    end
end

%% generate infeasible (NOT SPARSE,SLOW AS A RESULT)
% A'y = 0, y \in K*, b'*y = -1
if (gen_infeasible)
    z = randn(m,1);
    y = proj_dual_cone(z,K); % y = s - z;
    A = randn(m,n);
    
    A = A - ((A'*y)*y'/norm(y)^2)'; % dense...
    
    b = randn(m,1);
    b = -b / (b'*y);
    
    data.A = sparse(A);
    data.b = b;
    data.c = randn(n,1);
    
    params.scale = 0.5;
    %indirect
    if(run_indirect) [xi,yi,si,infoi] = scs_indirect(data,K,params); end
    % direct:
    if (run_direct) [xd,yd,sd,infod] = scs_direct(data,K,params); end
    
    % cvx:
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
    
    if(run_scs_matlab)
        params.use_indirect = true;
        [xi_m,yi_m,si_m,infoi_m] = scs_matlab(data,K,params);
        params.use_indirect = false;
        [xd_m,yd_m,sd_m,infod_m] = scs_matlab(data,K,params);
    end
    
end
%% generate unbounded (NOT SPARSE,SLOW AS A RESULT)
% Ax + s = 0, s \in K, c'*x = -1
if(gen_unbounded)
    z = randn(m,1);
    s = proj_cone(z,K);
    A = randn(m,n);
    x = randn(n,1);
    A = A - (s + A*x)*x'/(norm(x)^2); % dense...
    c = randn(n,1);
    c = - c / (c'*x);
    
    data.A = sparse(A);
    data.b = randn(m,1);
    data.c = c;
    
    params.scale = 0.5;
    %indirect
    if(run_indirect) [xi,yi,si,infoi] = scs_indirect(data,K,params); end
    % direct:
    if (run_direct) [xd,yd,sd,infod] = scs_direct(data,K,params); end
    
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
    
    if(run_scs_matlab)
        params.use_indirect = true;
        [xi_m,yi_m,si_m,infoi_m] = scs_matlab(data,K,params);
        params.use_indirect = false;
        [xd_m,yd_m,sd_m,infod_m] = scs_matlab(data,K,params);
    end
end
