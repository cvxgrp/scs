clear all

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run.')
disp('It may also crash/run out of memory.')
disp('Set run_sdpt3 = false if you just want to run scs.')
disp('------------------------------------------------------------')

run ../../matlab/cvx_install_scs.m

save_results = false;
run_sdpt3 = false;
run_scs = true;

ns = [100,500,1000];
ms = ns; % square matrices, but doesn't have to be

time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';

for i = 1:length(ns)
    seedstr = sprintf('scs_rpca_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))
    
    n = ns(i);
    m = ms(i);
    r = 10; % rank
    
    L1 = randn(m,r);
    L2 = randn(r,n);
    L = L1*L2;
    S = 10*sprandn(m,n,0.1);
    M = L + S;
    kap = sum(norms(S,1));
    
    %%
    if run_scs
        
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('SCALE',5)
        variables Lc(m,n) Sc(m,n)
        dual variable Yc
        minimize(norm_nuc(Lc))
        sum(norms(Sc,1)) <= kap
        Yc:Lc + Sc == M;
        output = evalc('cvx_end')
        toc
        
        scs_direct.L{i} = Lc;
        scs_direct.obj(i) = norm_nuc(Lc);
        scs_direct.output{i} = output;
        
        
        if (save_results); save('data/rpca_scs_direct', 'scs_direct'); end
        
        %%
        
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('USE_INDIRECT',1,'SCALE',5)
        variables Lc(m,n) Sc(m,n)
        dual variable Yc
        minimize(norm_nuc(Lc))
        sum(norms(Sc,1)) <= kap
        Yc:Lc + Sc == M;
        output = evalc('cvx_end')
        toc
        
        scs_indirect.L{i} = Lc;
        scs_indirect.obj(i) = norm_nuc(Lc);
        scs_indirect.output{i} = output;
        
        if (save_results); save('data/rpca_scs_indirect', 'scs_indirect'); end
        
    end
    %%
    if run_sdpt3
        try
            tic
            cvx_begin
            cvx_solver sdpt3
            variables Lt(m,n) St(m,n)
            dual variable Yt
            minimize(norm_nuc(Lt))
            sum(norms(St,1)) <= kap
            Yt:Lt + St == M;
            output = evalc('cvx_end')
            toc
            
            cvx.L{i} = Lt;
            cvx.obj(i) = norm_nuc(Lt);
            timing = regexp(output, time_pat_cvx, 'names');
            cvx.time{i} = str2num(timing.total);
            cvx.output{i} = output;
            cvx.err{i} = 0;
            
        catch err
            cvx.time{i} = toc;
            cvx.err{i} = err;
        end
        
        if (save_results); save('data/rpca_cvx', 'cvx'); end
        
    end
end
