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

ns = [5000, 50000, 100000];
ms = [50, 500, 1000];

density = 0.01;

time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';

for i = 1:length(ns)
    seedstr = sprintf('scs_portfolio_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))
    
    n = ns(i);
    m = ms(i);
    
    mu = exp(randn(n,1));
    D = sqrt(2*rand(n,1));
    F = sprandn(n,m,density);
    gamma = 10;
    
    portfolio_prob.F{i} = F;
    portfolio_prob.D{i} = D;
    portfolio_prob.mu{i} = mu;
    portfolio_prob.n{i} = n;
    portfolio_prob.m{i} = m;
    portfolio_prob.gamma{i} = gamma;
    
    if (save_results); save('data/portfolio_prob', 'portfolio_prob','-v7.3'); end
    
    %%
    if run_scs
        
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver scs_matlab'
        variable x(n)
        maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
        sum(x) == 1
        x >= 0
        output = evalc('cvx_end')
        toc
        
        scs_direct.x{i} = x;
        scs_direct.x_viol{i} = min(x);
        scs_direct.budget_viol{i} = abs(1-sum(x));
        scs_direct.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
        scs_direct.output{i} = output;
        
        if (save_results); save('data/portfolio_scs_direct', 'scs_direct'); end
        
        %%
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',2)
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver scs_matlab'
        variable x(n)
        maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
        sum(x) == 1
        x >= 0
        output = evalc('cvx_end')
        toc
        
        scs_indirect.x{i} = x;
        scs_indirect.x_viol{i} = min(x);
        scs_indirect.budget_viol{i} = abs(1-sum(x));
        scs_indirect.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
        scs_indirect.output{i} = output;
        
        
        if (save_results); save('data/portfolio_scs_indirect', 'scs_indirect'); end
        
    end
    %%
    if run_sdpt3
        try
            tic
            cvx_begin
            cvx_solver sdpt3
            variable x(n)
            maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
            sum(x) == 1
            x >= 0
            output = evalc('cvx_end')
            toc
            
            cvx.x{i} = x;
            cvx.x_viol{i} = min(x);
            cvx.budget_viol{i} = abs(1-sum(x));
            cvx.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
            timing = regexp(output, time_pat_cvx, 'names');
            cvx.time{i} = str2num(timing.total);
            cvx.output{i} = output;
            cvx.err{i} = 0;
            
        catch err
            cvx.time{i} = toc;
            cvx.err{i} = err;
        end
        
        if (save_results); save('data/portfolio_cvx', 'cvx'); end
        
    end
end