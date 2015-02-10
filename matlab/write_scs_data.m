function write_scs_data(data,K,params,name)
% this function takes in problem data and writes out an
% uncompressed file that can be read by demo_direct and
% demo_indirect (compiled scs binary files)
% this is just for testing/debugging as the files
% for some problems can be prohibitively large

% set default params if not present:
max_iters = 2500; % maximum num iterations for admm
eps = 1e-3; % quitting tolerances
alpha = 1.8; % relaxation parameter (alpha = 1 is unrelaxed)
normalize = 1; % 1 = true
verbose = 1; % 1 = true
rho_x = 1e-3;
scale = 5;
cg_rate = 2;
if ~isfield(params,'max_iters');params.max_iters = max_iters;end
if ~isfield(params,'eps');params.eps = eps;end
if ~isfield(params,'alpha');params.alpha = alpha;end
if ~isfield(params,'normalize');params.normalize = normalize;end
if ~isfield(params,'verbose');params.verbose = verbose;end
if ~isfield(params,'rho_x');params.rho_x = rho_x;end
if ~isfield(params,'scale');params.scale = scale;end
if ~isfield(params,'cg_rate');params.cg_rate = cg_rate;end

%% set cone data if not present:
if ~isfield(K,'f');K.f = 0;end
if ~isfield(K,'l');K.l = 0;end
if ~isfield(K,'q');K.q = [];end
if ~isfield(K,'s');K.s = [];end
if ~isfield(K,'ep');K.ep = 0;end
if ~isfield(K,'ed');K.ed = 0;end
if ~isfield(K,'p');K.p = [];end

n = length(data.c);
m = size(data.A,1);
data.A = sparse(data.A);

% col-compressed A
[i,~,s] = find(sparse(data.A));
i = i-1;
tmp = full(sum(data.A~=0));
pw = [0 cumsum(tmp)];

%Q=sparse([zeros(n) data.A' data.c;
%    -data.A zeros(m,m) data.b;
%    -data.c' -data.b' 0]);
%W=sparse([speye(n+m+1) Q';Q -speye(n+m+1)]);

delete(name);
fi = fopen(name,'w');
fprintf(fi,'%u ',n);fprintf(fi,'%u ',m);
fprintf(fi,'\n');
fprintf(fi,'%u ',K.f);fprintf(fi,'%u ',K.l);fprintf(fi,'%u ',length(K.q));fprintf(fi,'%u ',length(K.s));fprintf(fi,'%u ',K.ep);fprintf(fi,'%u ',K.ed);fprintf(fi,'%u ',length(K.p));
fprintf(fi,'\n');
fprintf(fi,'%u ',params.max_iters);
fprintf(fi,'\n');
fprintf(fi,'%u ',params.verbose); fprintf(fi,'%u ',params.normalize);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.alpha);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.eps);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.rho_x);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.scale);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.cg_rate);
fprintf(fi,'\n');
fprintf(fi,'%u ',K.q');
fprintf(fi,'\n');
fprintf(fi,'%u ',K.s');
fprintf(fi,'\n');
fprintf(fi,'%u ',K.p');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',full(data.b)');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',full(data.c)');
fprintf(fi,'\n');

fprintf(fi,'%u ',pw);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');

fclose(fi);
end
