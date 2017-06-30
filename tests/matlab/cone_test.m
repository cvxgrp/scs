clc
clear all;

K.q = 4;

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

params.eps          = 1e-8;             
params.nominal      = 0;
params.do_super_scs = 1;
params.alpha        = 1.5;
params.scale        = 1;
params.verbose      = 1;
params.normalize    = 1;
params.direction    = 100;
params.beta         = 0.5;
params.c1           = 1.0 - 1e-4;
params.c_bl         = 0.999;
params.k0           = 0;
params.k1           = 1;
params.k2           = 1;
params.ls           = 10;
params.sigma        = 1e-2;
params.thetabar     = 0.1;
params.memory       = 10;
params.sse          = 0.999;
params.tRule        = 1;
params.do_record_progress = 1;
params.max_iters    = 1000;
params.rho_x        = 0.2;
[x2, y2, s2, info2] = superscsCversion(data, K, params);
[x1, y1, s1, info1] = scs_direct(data, K, params);




fprintf('|errx| = %g, |erry| = %g, |errs| = %g\n', ...
    norm(x1 - x2, Inf), norm(y1 - y2, Inf), norm(s1 - s2, Inf));
assert(norm(x1 - x2, Inf)<1e-7,'x');
if (all(~isnan(y1)) && all(~isnan(y2))), assert(norm(y1 - y2, Inf)<1e-6,'y'); end
assert(norm(s1 - s2, Inf)<1e-7,'z');

% info1.iter - info2.iter
% assert((info1.resPri-info2.resPri)/info2.resPri<1e-1, 'resPri');
% assert((info1.relGap-info2.relGap)/info2.relGap<1e-1, 'relative gap');
% info1.resDual - info2.resDual

if strcmp('Solved', info1.status)==1,
    tol = max([info1.resPri,info1.resDual,info1.relGap]);
    assert(tol < params.eps, 'inaccurate solution')
end
%%

A(1,1) = 0.3; A(4,1) = -0.5;
A(2,2) = 0.7; A(4,2) = 0.9; A(3,3) = 0.2;
A = sparse(A);

b = [0.2; 0.1; -0.1; 0.1];
c = [1;-2;-3];

n = size(A,2);
cvx_begin
cvx_solver scs
cvx_solver_settings('eps', 1e-8, 'do_super_scs', 1, 'rho_x', 1, 'direction', 100, 'memory', 50 );
variable x(n);
dual variable y;
minimize( x'*x + c' * x );
subject to
y : A * x <= b;
cvx_end