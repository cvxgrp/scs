function run_pnorm_ex(params)

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run.')
disp('It may also crash/run out of memory.')
disp('Set run_cvx = false if you just want to run scs.')
disp('------------------------------------------------------------')

save_results = false;
run_cvx = true;
cvx_use_solver = 'sdpt3';
run_scs_direct = true;
run_scs_indirect = true;

if nargin==1
    if isfield(params,'save_results');  save_results = params.save_results;end
    if isfield(params,'run_cvx');       run_cvx = params.run_cvx;end
    if isfield(params,'cvx_use_solver');cvx_use_solver = params.cvx_use_solver;end
    if isfield(params,'run_scs_direct');run_scs_direct = params.run_scs_direct;end
    if isfield(params,'run_scs_indirect');run_scs_indirect = params.run_scs_indirect;end
end

ns = [1000, 10000, 100000];
ms = ceil(ns/2);

density = 0.1;

for i = 1:length(ns)
    seedstr = sprintf('scs_pnorm_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))
    
    n = ns(i);
    m = ms(i);
    pow = pi;
    G = sprandn(m,n,density);
    f = randn(m,1) * n * density;
    
    n=ns(i);
    m=ms(i);
    
    %%
    if (run_scs_direct || run_scs_indirect)
        %% scs, with power cone formulation
        Anz = nnz(G) + n + 1 + 3*n;
        data.A = sparse([],[],[],m+3*n,2*n+1,Anz);
        data.A(1:m+1,:) = [G sparse(m,n) sparse(m,1); sparse(1,n) ones(1,n) -1];
        for j=1:n
            ej = sparse(n,1);ej(j) = -1;
            data.A(m+1+(j-1)*3+1:m+1+j*3,:) = [sparse(1,n) ej' 0; sparse(1,2*n) -1;  ej' sparse(1,n) 0];
        end
        data.c = [zeros(2*n,1); 1];
        data.b = [f; 0 ; zeros(3*n,1)];
        K = struct('f', m+1, 'p', ones(n,1) / pow);
        
        params.eps = 1e-6;
        params.scale = 1;
        params.cg_rate = 1.5;
        params.do_super_scs = 1;
        params.memory = 50;
        params.normalize = 1;
        params.direction = 100;
        params.rho_x = 1;
        if (run_scs_direct)
            if (save_results)
                [out,x_scs,y_scs,s_scs,info] = evalc('scs_direct(data, K, params)');
                out
                save('data/pnorm_scs_direct', 'out');
            else
                [x_scs,y_scs,s_scs,info] = scs_direct(data, K, params);
            end
        end
        
        if (run_scs_indirect)
            if (save_results)
                [out,x_scs,y_scs,s_scs,info] = evalc('scs_indirect(data, K, params)');
                out
                save('data/pnorm_scs_indirect', 'out');
            else
                [x_scs,y_scs,s_scs,info] = scs_indirect(data, K, params);
            end
        end
        
    end
    %%
    if run_cvx
        try
            tic
            cvx_begin
            cvx_solver scs
            cvx_solver_settings('eps',1e-4,'verbose',2,'do_super_scs', 1, 'rho_x', 1, 'direction', 100, 'memory', 40 );
    
            variable x(n)
            minimize(norm(x, pow))
            subject to
            G*x == f
            if (save_results)
                cvx.output{i} = evalc('cvx_end')
            else
                cvx_end
            end
            toc
            
        catch err
            err
            cvx.err{i} = err;
        end
        
        if (save_results); save(sprintf('data/pnorm_cvx_%s',cvx_use_solver), 'cvx'); end;
        
    end
end
