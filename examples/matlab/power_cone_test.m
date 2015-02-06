clear all;
randn('seed',0);
n = 50;
pow = pi/2;
const = 1.5;
a = randn(n,1);

%% cvx, using SOCP formulation
cvx_begin
%cvx_solver('sedumi')
variable x(n)
minimize(norm(x, pow))
subject to
a'*x == const
cvx_end

%% scs, with power cone formulation
data.A = [a' zeros(1,n) 0; zeros(1,n) ones(1,n) -1];
for i=1:n
    ei = zeros(n,1);ei(i) = -1;
    data.A = [data.A; zeros(1,n) ei' 0; zeros(1,2*n) -1;  ei' zeros(1,n) 0];
end
data.A = sparse(data.A);
data.c = [zeros(2*n,1); 1];
data.b = [const; 0 ; zeros(3*n,1)];
K = struct('f', 2, 'p', ones(n,1) / pow);

params.eps = 1e-6;
[x_scs,y_scs,s_scs,info] = scs_indirect(data, K, params);

%params.gen_plots = true;
%params.use_indirect = true;
%[x_scs,y_scs,s_scs,info] = scs_matlab(data, K, params);

norm(x-x_scs(1:n))