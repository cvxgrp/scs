write_scs_data(data,K,params,name)
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
if ~isfield(params,'max_iters');params.max_iters = max_iters;end
if ~isfield(params,'eps');params.eps = eps;end
if ~isfield(params,'alpha');params.alpha = alpha;end
if ~isfield(params,'normalize');params.normalize = normalize;end
if ~isfield(params,'verbose');params.VERBOSE = VERBOSE;end

%% set cone data if not present:
if ~isfield(K,'f');K.f = 0;end
if ~isfield(K,'l');K.l = 0;end
if ~isfield(K,'q');K.q = [];end
if ~isfield(K,'s');K.s = [];end
if ~isfield(K,'ep');K.ep = 0;end
if ~isfield(K,'ed');K.ed = 0;end

n = length(data.c);
m = size(data.A,1);
data.A = sparse(data.A);

%{
% symmetrize SD cone matrices:
[mm,nn]=size(data.A);
idx = K.f + K.l + sum(K.q);
for i=1:size(K.s)
    for j=1:nn
        work = data.A(idx+1:idx+K.s(i)^2, j);
        work = sparse(reshape(work,K.s(i),K.s(i)));
        if any(any(work~=work'))
            %fprintf('warning: symmetrizing A\n')
            work = (work+work')/2;
            data.A(idx+1:idx+K.s(i)^2, j) = sparse(work(:));
        end
    end
    
    work = data.b(idx+1:idx+K.s(i)^2);
    work = sparse(reshape(work,K.s(i),K.s(i)));
    if any(any(work~=work'))
        %fprintf('warning: symmetrizing b\n')
        work = (work+work')/2;
        data.b(idx+1:idx+K.s(i)^2) = sparse(work(:));
    end
    
    idx = idx + K.s(i)^2;
end
clear work;
%}
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
fprintf(fi,'%u ',K.f);fprintf(fi,'%u ',K.l);fprintf(fi,'%u ',length(K.q));fprintf(fi,'%u ',length(K.s));fprintf(fi,'%u ',K.ep);fprintf(fi,'%u ',K.ed);
fprintf(fi,'\n');
fprintf(fi,'%u ',params.max_iters);
fprintf(fi,'\n');
fprintf(fi,'%u ',params.verbose); fprintf(fi,'%u ',params.normalize);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.alpha);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.eps);
fprintf(fi,'\n');
fprintf(fi,'%u ',K.q');
fprintf(fi,'\n');
fprintf(fi,'%u ',K.s');
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

% triplet A
%{
[i,j,s] = find(data.A);
i = i-1;j=j-1;
NNZ=length(i);
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',j');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');
%}

fclose(fi);
