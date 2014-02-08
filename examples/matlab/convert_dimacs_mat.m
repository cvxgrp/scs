close all; clear all
input_path = 'DIMACS_mat_files'
output_path = 'DIMACS'

mkdir(output_path)
cd(input_path)

tests = dir('*.mat');
for i=1:length(tests)
    szes(i) = tests(i).bytes;
end
[szes, solve_order] = sort(szes);
params = struct('VERBOSE', 1, 'EPS_ABS', 1e-3, 'MAX_ITERS', 2500, 'CG_MAX_ITS', 15);
%%
N = length(tests);
for ii = 1:N
    i = solve_order(ii); %% solve in increasing order of size
    
    clear A At b c K
    test_name = tests(i).name;
    f = [input_path '/' test_name];
    test_name = test_name(1:end-4); % strip .mat
    fprintf('converting test %i out of %i : %s\n', ii, N, test_name);
    
    load(f)
    
    [m1,n1] = size(b);
    [m2,n2] = size(c);
    if m1 == 1,
        b = b';
    end
    if m2 == 1,
        c = c';
    end
    if (exist('A','var'))
        data = struct('A', sparse(A'), 'b', full(c), 'c', -full(b));
    else
        data = struct('A', sparse(At), 'b', full(c), 'c', -full(b));
    end
    ccones = {'f','l','q','s'};
    for j = 1:length(ccones)
        if ~isfield(K,ccones{j})
            K.(ccones{j}) = [];
        end
    end
    if isempty(K.f)
        K.f = 0;
    end
    if isempty(K.l)
        K.l = 0;
    end
    if (K.q == 0)
        K.q = [];
    end
    if (K.s == 0)
        K.s = [];
    end
    
    cone = struct('f', 0, 'l', K.l, 'q', K.q', 's', K.s');
    
    [m1,n1] = size(cone.q);
    if m1 == 1,
        cone.q = cone.q';
    end
    [m1,n1] = size(cone.s);
    if m1 == 1,
        cone.s = cone.s';
    end
    if (data.b == 0); data.b = zeros(size(data.A,1),1); end
    if (data.c == 0); data.c = zeros(size(data.A,2),1); end
    
    if (isfield(K,'r') && K.r ~= 0)
        disp('rotated lorentz cones not currently supported');
    else
        cd ..
        write_coneOS_data_sparse(data,cone,params,sprintf('%s/%s',output_path,test_name));
        cd(input_path)
    end
end
cd ..
