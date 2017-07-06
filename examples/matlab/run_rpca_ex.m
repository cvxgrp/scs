function run_rpca_ex(params)

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run.')
disp('It may also crash/run out of memory.')
disp('Set run_cvx = false if you just want to run scs.')
disp('------------------------------------------------------------')

save_results = false;
run_cvx = true;
cvx_use_solver = 'scs';
run_scs_direct = true;
run_scs_indirect = true;

if nargin==1
    if isfield(params,'save_results');  save_results = params.save_results;end
    if isfield(params,'run_cvx');       run_cvx = params.run_cvx;end
    if isfield(params,'cvx_use_solver');cvx_use_solver = params.cvx_use_solver;end
    if isfield(params,'run_scs_direct');run_scs_direct = params.run_scs_direct;end
    if isfield(params,'run_scs_indirect');run_scs_indirect = params.run_scs_indirect;end
end

ns = [100, 500, 1000];
ms = ns; % square matrices, but doesn't have to be
density = 0.1;

for i = 1:length(ns)
    seedstr = sprintf('scs_rpca_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))
    
    n = ns(i);
    m = ms(i);
    r = 10; % rank
    
    L1 = randn(m,r);
    L2 = randn(r,n);
    L = L1*L2;
    S = 10*sprandn(m,n,density);
    M = L + S;
    kap = sum(sum(abs(S)));
    
    %%
    if run_scs_direct
        
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('eps',1e-3,'scale',1)
        variables Lc(m,n) Sc(m,n)
        dual variable Yc
        minimize(norm_nuc(Lc))
        sum(norms(Sc,1)) <= kap
        Yc:Lc + Sc == M;
        if (save_results)
            output = evalc('cvx_end')
        else
            output='';
            cvx_end
        end
        toc
        
        scs_direct.L{i} = Lc;
        scs_direct.obj(i) = norm_nuc(Lc);
        scs_direct.output{i} = output;
        
        
        if (save_results); save('data/rpca_scs_direct', 'scs_direct'); end
    end
    %%
    if run_scs_indirect
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('use_indirect',1,'eps',1e-3,'scale',1,'cg_rate',1.5)
        variables Lc(m,n) Sc(m,n)
        dual variable Yc
        minimize(norm_nuc(Lc))
        sum(norms(Sc,1)) <= kap
        Yc:Lc + Sc == M;
        if (save_results)
            output = evalc('cvx_end')
        else
            output='';
            cvx_end
        end
        toc
        
        scs_indirect.L{i} = Lc;
        scs_indirect.obj(i) = norm_nuc(Lc);
        scs_indirect.output{i} = output;
        
        if (save_results); save('data/rpca_scs_indirect', 'scs_indirect'); end
    end
    %%
    if run_cvx
        try
            tic
            cvx_begin
            cvx_solver(cvx_use_solver)
            variables Lt(m,n) St(m,n)
            dual variable Yt
            minimize(norm_nuc(Lt))
            sum(norms(St,1)) <= kap
            Yt:Lt + St == M;
            if (save_results)
                output = evalc('cvx_end')
            else
                output='';
                cvx_end
            end
            toc
            
            cvx.L{i} = Lt;
            cvx.obj(i) = norm_nuc(Lt);
            cvx.output{i} = output;
            cvx.err{i} = 0;
            
        catch err
            err
            cvx.err{i} = err;
        end
        
        if (save_results); save(sprintf('data/rpca_cvx_%s',cvx_use_solver), 'cvx'); end
        
    end
end
