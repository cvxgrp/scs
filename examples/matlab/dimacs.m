close all; clear all
cd 'DIMACS_mat_files'

%run ../../coneOSsparse/matlab/install_scs_cvx.m
%copyfile('../../coneOSsparse/matlab/scs_direct.m*','.');

cvx_on = true;
%cvx_use_solver = 'sdpt3';

scs_on = false;
save_results = false;

tests = dir('*.mat');
params = struct('VERBOSE', 1, 'EPS_ABS', 1e-5, 'MAX_ITERS', 10000);
for i=1:length(tests)
    szes(i) = tests(i).bytes;
end
[szes, solve_order] = sort(szes);


%time_pat_scs = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
%iter_pat_scs = {'(?<iter>[\d]+)\|'};

%%

if exist('../data/dimacs_sdpt3.mat')
    load('../data/dimacs_sdpt3.mat')
end

N = length(tests);

[ ver, isoctave, fs, ps ] = cvx_version;
addpath(strcat(fs,'/sdpt3/'));
addpath(strcat(fs,'/sdpt3/Solver'));
addpath(strcat(fs,'/sdpt3/Solver/Mexfun'));
addpath(strcat(fs,'/sdpt3/HSDSolver'));


for ii = 1:N
    i = solve_order(ii); %% solve in increasing order of size
    
    %clear A At b c K
    clear blk A C b output obj X y Z info runhist
    test_name = tests(i).name;
    f = ['DIMACS_mat_files/' test_name];
    test_name = test_name(1:end-4); % strip .mat
    fprintf('running test %i out of %i : %s\n', ii, N, test_name);
    
    %load(f)
    %{
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
    %}
    
    if cvx_on
        
        if (~exist('dimacs_sdpt3') || ~isfield(dimacs_sdpt3,test_name))
            try
                
                %{
                [m,n] = size(data.A);
                
                cvx_begin %quiet
                cvx_solver(cvx_use_solver)
                %cvx_solver_settings('MAX_ITERS',2500,'EPS',1e-3)
                variables xcvx(n) scvx(m)
                dual variable zcvx
                minimize(data.c'*xcvx)
                zcvx: data.A*xcvx + scvx == data.b
                scvx(1:cone.f) == 0
                scvx(cone.f+1:cone.f+cone.l) >= 0
                idx=cone.f+cone.l;
                if isempty(idx) idx = 0; end
                for kk =1:length(cone.q)
                    norm(scvx(idx + 2: idx + cone.q(kk))) <= scvx(idx + 1)
                    idx = idx + cone.q(kk);
                end
                for kk =1:length(cone.s)
                    reshape(scvx(idx+1:idx + cone.s(kk)^2),cone.s(kk),cone.s(kk)) == semidefinite(cone.s(kk));
                    idx = idx + cone.s(kk)^2;
                end
                output = evalc('cvx_end')
                
                sdpt3.(test_name).obj = cvx_optval;
                timing = regexp(output, time_pat_cvx, 'names');
                sdpt3.(test_name).time = str2num(timing.total);
                sdpt3.(test_name).output = output;
                %}
                
                [blk,A,C,b] = read_sedumi(f);
                %[output,obj,X,y,Z,info,runhist] = evalc('sdpt3(blk,A,C,b)');
                [obj,X,y,Z,info,runhist] = sdpt3(blk,A,C,b);

                output
                
                timing = regexp(output, time_pat_cvx, 'names');
                dimacs_sdpt3.(test_name).time = str2num(timing.total);
                dimacs_sdpt3.(test_name).info = info;
                dimacs_sdpt3.(test_name).output = output;
                
            catch err
                dimacs_sdpt3.(test_name).err = err;
            end
            if (save_results); save('../data/dimacs_sdpt3', 'dimacs_sdpt3','-v7.3'); end;
            
        end
        
    end
    
    if scs_on
        if (isfield(K,'r') && K.r ~= 0)
            scs_error(test_name) = 'rotated lorentz cones not currently supported';
            scs_x = nan;
            scs_objval = nan;
        else
            
            tic
            [output, x_m, z_m] = evalc('scs_direct(data,cone,params);');
            output
            data.b'*z_m %% because scs targets dual of DIMACs formulations
            toc
            
            scs_NNZA.(test_name) = nnz(data.A);
            scs_x.(test_name) = x_m;
            scs_z.(test_name) = z_m;
            timing = regexp(output, time_pat_scs, 'names');
            scs_times.(test_name) = str2num(timing.total);
            tmp = regexp(output,iter_pat_scs, 'names');
            scs_iters.(test_name) = str2num(tmp{1}(end).iter) + 1;
            scs_dobjval.(test_name) = data.b'*z_m;
            scs_pobjval.(test_name) = data.c'*x_m;
            scs_output.(test_name) = output;
            
            if save_results; save('../data/dimacs_scs','scs','-v7.3'); end
        end
    end
    
end
delete 'scs_direct.m*'
cd ..