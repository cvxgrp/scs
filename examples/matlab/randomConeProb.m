clear all;close all;
addpath('../../matlab')
cd '../../matlab'; make_scs; cd '../examples/matlab';

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

% set cone sizes (ep = ed = 0 if you want to compare against cvx):
K = struct('f',100,'l',150,'q',[2;3;4;5;6;7;8;9;10;5;6;100;0;1],'s',[0,1,2,3,4,5],'ep',5,'ed',5)
density = 0.01; % A matrix density

m = getConeDims(K);
n = round(m/3);
params = struct('EPS',1e-3, 'NORMALIZE',1,'SCALE',5,'CG_RATE',2);

%% generate primal-dual feasible cone prob:
% Ax + s = b, s \in K, A'y + c = 0, y \in K*, s'y = 0
if (gen_feasible)
    z = randn(m,1);
    z = symmetrizeSDP(z,K); % for SD cones
    y = proj_dual_cone(z,K); % y = s - z;
    s = y - z; %s = proj_cone(z,K);
    
    
    A = sprandn(m,n,density);
    x = randn(n,1);
    c = -A'*y;
    b = A*x + s;
    nnz(A)
    
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
end

%% generate infeasible (NOT SPARSE,SLOW AS A RESULT)
% A'y = 0, y \in K*, b'*y = -1
if (gen_infeasible)
    z = randn(m,1);
    z = symmetrizeSDP(z,K); % for SD cones
    y = proj_dual_cone(z,K); % y = s - z;
    A = randn(m,n);
    
    A = A - ((A'*y)*y'/norm(y)^2)'; % dense...
    
    b = randn(m,1);
    b = -b / (b'*y);
    
    data.A = sparse(A);
    data.b = b;
    data.c = randn(n,1);
    
    params.SCALE = 0.5;
    %indirect
    if(run_indirect) [xi,yi,si,infoi] = scs_indirect(data,K,params); end
    % direct:
    if (run_direct) [xd,yd,sd,infod] = scs_direct(data,K,params); end
    
    % cvx:
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
    
end
%% generate unbounded (NOT SPARSE,SLOW AS A RESULT)
% Ax + s = 0, s \in K, c'*x = -1
if(gen_unbounded)
    z = randn(m,1);
    z = symmetrizeSDP(z,K); % for SD cones
    s = proj_cone(z,K);
    A = randn(m,n);
    x = randn(n,1);
    A = A - (s + A*x)*x'/(norm(x)^2); % dense...
    c = randn(n,1);
    c = - c / (c'*x);
    
    data.A = sparse(A);
    data.b = randn(m,1);
    data.c = c;
    
    params.SCALE = 0.5;
    %indirect
    if(run_indirect) [xi,yi,si,infoi] = scs_indirect(data,K,params); end
    % direct:
    if (run_direct) [xd,yd,sd,infod] = scs_direct(data,K,params); end
    
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
end
