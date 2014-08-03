clear all
make_scs

fs = cvx_where;
if (strcmp(fs(end),filesep));fs = fs(1:end-1);end

scs_install_path = strcat(fs,filesep,'scs');
mkdir(scs_install_path);

copyfile(strcat('scs_direct.',mexext),scs_install_path);
copyfile(strcat('scs_indirect.',mexext),scs_install_path);

%copyfile('cvx_scs.m', strcat(fs,filesep,'shims',filesep,'cvx_scs.m'));
copyfile('scs.m', scs_install_path);
cvx_setup

%{
n = 10;m = 20;
A = randn(m,n);
b = randn(m,1);
cvx_begin
cvx_solver 'scs'
%cvx_solver_settings('MAX_ITERS',10)
%cvx_solver_settings('SCALE',5)
variable x(n)
minimize(norm(A*x - b))
cvx_end
%}