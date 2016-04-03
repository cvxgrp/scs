function [x,y,s,info] = scs_matlab(data,K,params)
% cone solver, solves:
%
% min. c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% K is product of cones in this exact order:
%
%       zero cone (dual to free cone)
%       linear cone
%       second-order cones
%       semidefinite cones
%       exponential cones
%       dual exponential cones
%
% ** The rows of A must be in the same order as the cones above **
%
% data must consist of data.A, data.b, data.c, where A,b,c used as above
%
% cone struct must consist of:
%
%       cone.f, length of zero cone (for equality constraints)
%       cone.l, length of lp cone
%       cone.q, array of SOC lengths
%       cone.s, array of SD lengths
%       cone.ep, num of exp cones
%       cone.ed, num of dual exp cones
%       cone.p, array of primal + dual power cone params
%
% cone struct is only used in proj_cone, to add other cones
% simply add the relevant size data to the cone struct and edit the
% proj_cone method to include projection onto the new cone

% params struct consists of the following fields
% (here set to default settings):
gen_plots = false;  % generate convergence plots
max_iters = 2500;   % maximum num iterations for admm
eps = 1e-3;         % quitting tolerances
alpha = 1.5;        % relaxation parameter (alpha = 1 is unrelaxed)
normalize = 1;      % heuristic normalization procedure
scale = 1;          % heuristic re-scaline procedure
rho_x = 1e-3;       % x equality rescaling

% conjugate gradient (CG) settings:
use_indirect = false;   % use conjugate gradient rather than direct method
cg_rate = 2;            % max CG quitting tolerance
extra_verbose = false;  % CG prints summary

% constants
UNDET_TOL = 1e-9;   % tol for undetermined solution (tau = kappa = 0)
PRINT_INTERVAL = 100;
CONVERGED_INTERVAL = 20;
%%
if nargin==3
    if isfield(params,'gen_plots');gen_plots = params.gen_plots;end
    if isfield(params,'max_iters');max_iters = params.max_iters;end
    if isfield(params,'eps');eps = params.eps;end
    if isfield(params,'alpha');alpha = params.alpha;end
    if isfield(params,'normalize');normalize = params.normalize;end
    if isfield(params,'scale');scale = params.scale;end
    if isfield(params,'rho_x');rho_x = params.rho_x;end
    % CG:
    if isfield(params,'use_indirect');use_indirect = params.use_indirect;end
    if isfield(params,'cg_rate');cg_rate = params.cg_rate;end
    if isfield(params,'extra_verbose');extra_verbose = params.extra_verbose;end
end
%%
% Q matrix (as in paper):
%{
Q = sparse([zeros(n) data.A' data.c;
           -data.A zeros(m,m) data.b;
           -data.c' -data.b' 0]);
%}

K = validateCone(K);

n = length(data.c);
m = length(data.b);
l=n+m+1;

nm_b = norm(data.b);
nm_c = norm(data.c);

work = struct();
if (normalize)
    [data, work] = normalize_data(data, K, scale, work); % have to do this since matlab pass by val
    D = work.D;
    E = work.E;
    sc_b = work.sc_b;
    sc_c = work.sc_c;
else
    scale = 1;
    D = ones(m,1);
    E = ones(n,1);
    sc_b = 1;
    sc_c = 1;
end


%%
if use_indirect
    work.M = 1 ./ diag(rho_x*speye(n) + data.A'*data.A); % pre-conditioner
else
    W=sparse([rho_x*speye(n) data.A';data.A -speye(m)]);
    disp('Factorization')
    try
        work.P=amd(W);
        [work.L,work.d] = ldlsparse(W,work.P);
        work.L=work.L+speye(n+m);
    catch ldlerror
        disp('WARNING: LDLSPARSE ERROR, using MATLAB LDL instead (this is slower).')
        [work.L,work.d,work.P] = ldl(W,'vector');
    end
end

h = [data.c;data.b];
[g, gTh, ~] = solve_for_g(work,data,h,n,m,rho_x,use_indirect,cg_rate,extra_verbose);

% u = [x;z;tau], v = [y;s;kappa]
fprintf('Iter:\t      pres      dres       gap      pobj      dobj   unb_res   inf_res   kap/tau  time (s)\n');

% for warm-starting:
if (nargin==3 && isfield(params,'warm_xy'))
    u = [params.warm_xy;1];
    v = [zeros(n,1);data.b*u(end) - data.A*u(1:n);0];
    if (normalize)
        u(1:n) = u(1:n) .* (E * sc_b);
        u(n+1:n+m) = u(n+1:n+m) .* (D * sc_c);
        v(n+1:n+m) = v(n+1:n+m) ./ (D / (sc_b * scale));
    end
else
    u = zeros(l,1);u(end) = sqrt(l);
    v = zeros(l,1);v(end) = sqrt(l);
end

u_bar = u;
ut_bar = u;
v_bar = v;

tic
for i=0:max_iters-1
    % solve linear system
    u_prev = u;
    [ut, itn] = project_lin_sys(work, data, n, m, u, v, rho_x, i, use_indirect, cg_rate, extra_verbose,  h, g, gTh);
    %% K proj:
    rel_ut = alpha*ut+(1-alpha)*u;
    rel_ut(1:n) = ut(1:n); % don't relax 'x' variable
    u = rel_ut - v;
    u(n+1:n+m) = proj_dual_cone(u(n+1:n+m),K);
    u(l) = max(u(l),0);
    
    %% dual update:
    v = v + (u - rel_ut);
    
    % ergodic behavior
    u_bar = (u + u_bar * i) / (i+1);
    ut_bar = (ut + ut_bar * i) / (i+1);
    v_bar = (v + v_bar * i) / (i+1);

    %% convergence checking:
    tau = abs(u(end));
    kap = abs(v(end)) / (sc_b * sc_c * scale);
    
    x = u(1:n) / tau;
    y = u(n+1:n+m) / tau;
    s = v(n+1:n+m) / tau;
    
    err_pri = norm(D.*(data.A * x + s - data.b)) / (1 + nm_b) / (sc_b * scale);
    err_dual = norm(E.*(data.A' * y + data.c)) / (1 + nm_c) / (sc_c * scale);
    pobji = data.c' * x / (sc_c * sc_b * scale);
    dobji = -data.b' * y / (sc_c * sc_b * scale);
    
    gap = abs(pobji - dobji) / (1 + abs(pobji) + abs(dobji));
    
    if (data.c'*u(1:n) < 0)
        unb_res = norm(E.*data.c) * norm(D.*(data.A * u(1:n) + v(n+1:n+m))) / (-data.c'*u(1:n)) / scale;
    else
        unb_res = inf;
    end
    if (data.b'*u(n+1:n+m) < 0)
        inf_res = norm(D.*data.b) * norm(E.*(data.A' * u(n+1:n+m))) / (-data.b'*u(n+1:n+m)) / scale;
    else
        inf_res = inf;
    end
    idx = i+1;
    if use_indirect
        cg_its(idx) = itn;
        mults(idx) = 2+2*itn;
    end
    nms(idx,1) = err_pri;
    nms(idx,2) = err_dual;
    nms(idx,3) = gap;
    
    pathol(idx,1) = unb_res;
    pathol(idx,2) = inf_res;
    
    tau_i(idx) = ut(end);
    kap_i(idx) = v(end);
    pobj(idx) = pobji;
    dobj(idx) = dobji;
    
    sum_nm2(idx,1) = norm(u - u_prev);
    sum_nm2(idx,2) = norm(u - ut);
    
    solved = max(gap,max(err_pri,err_dual)) < eps;
    infeasible = inf_res < eps;
    unbounded = unb_res < eps;
    
    if (mod(i,CONVERGED_INTERVAL)==0 && (solved || infeasible || unbounded))
        break
    end
    
    if (mod(i,PRINT_INTERVAL)==0)
        ttime = toc;
        fprintf('%i:\t%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e\n',i,err_pri,err_dual,gap,pobji,dobji,unb_res,inf_res,kap/tau,ttime);
    end
end
if (i+1 == max_iters); i=i+1; end
ttime = toc;
fprintf('%i:\t%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e\n',i,err_pri,err_dual,gap,pobji,dobji,unb_res,inf_res,kap/tau,ttime);
if (solved)
    fprintf('\tc*x = %4f, -b*y = %4f, dist(s,K) = %.2e, dist(y,K*) = %.2e, s*y = %.2e\n', pobji,...
        dobji, norm(s - proj_cone(s,K)), norm(y - proj_dual_cone(y,K)), s'*y);
end
%%
tau = abs(u(end));
kap = abs(v(end)) / (sc_b * sc_c * scale);

x = u(1:n) / tau;
y = u(n+1:n+m) / tau;
s = v(n+1:n+m) / tau;

if (tau > UNDET_TOL && tau > kap) % this is different to Lieven
    status = 'Solved'
else
    x = nan(n,1);
    y = nan(m,1);
    s = nan(m,1);
    
    x_h = u(1:n);
    y_h = u(n+1:n+m);
    s_h = v(n+1:n+m);
    if norm((u+ut)/2)<=2*(UNDET_TOL*sqrt(l))
        status = 'Undetermined'
    elseif data.b'*y_h < data.c'*x_h
        status = 'Infeasible'
        y = -y_h * scale * sc_b * sc_c /(data.b'*y_h);
    else
        status = 'Unbounded'
        x = -x_h * scale * sc_b * sc_c /(data.c'*x_h);
        s = -s_h * scale * sc_b * sc_c /(data.c'*x_h);
    end
end
info.status = status;
info.iter = i;

info.resPri =err_pri;
info.resDual = err_dual;
info.relGap = gap;

if (normalize)
    y = y ./ (D * sc_c);
    x = x ./ (E * sc_b);
    s = s .* (D / (sc_b * scale));
    %[data, work] = unnormalize_data(data, scale, work);
end

%%
if use_indirect;
    fprintf('mean cg its: %4f\n', mean(cg_its));
end

if gen_plots
    figure();semilogy(nms(:,1));hold on;semilogy(nms(:,2),'r');semilogy(nms(:,3),'g');
    legend('pri resid','dual resid','gap');
    figure();plot(tau_i);hold on; plot(kap_i,'r');
    legend('tau','kappa')
    figure();plot(pobj);hold on;plot(dobj,'r');
    legend('primal obj','dual obj')
    figure();semilogy(sum_nm2(:,1));hold on;semilogy(sum_nm2(:,2),'r');
    legend('|u-uprev|','|u-ut|')
    figure();semilogy(pathol(:,1));hold on;semilogy(pathol(:,2),'r');
    legend('unb','inf');
    if use_indirect;
        figure();plot(cg_its);xlabel('k');ylabel('Conjugate Gradient Iterations');
        figure();semilogy(cumsum(mults),nms(:,1));hold on;semilogy(cumsum(mults),nms(:,2),'r');semilogy(cumsum(mults),nms(:,3),'g');
        legend('pri resid','dual resid','gap');
        xlabel('A multiplies');
    end
end
end
