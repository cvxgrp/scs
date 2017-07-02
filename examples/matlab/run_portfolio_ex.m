function run_portfolio_ex(params)

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run.')
disp('It may also crash/run out of memory.')
disp('Set run_cvx = false if you just want to run scs.')
disp('------------------------------------------------------------')

save_results = false;
run_cvx = false;
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

ns = [100000, 500000, 2500000];
ms = [100,    500,    2500];

density = 0.1;

for i = 1:length(ns)
    seedstr = sprintf('scs_portfolio_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))
    
    n = ns(i);
    m = ms(i);
    
    mu = exp(0.01*randn(n,1))-1; % returns
    D = rand(n,1)/10; % idiosyncratic risk
    F = sprandn(n,m,density)/10; % factor model
    gamma = 1;
    B = 1;
    %%
    if run_scs_direct
        
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('eps',1e-3,'scale',1,'do_super_scs',1,'direction',100,'memory',10,'k0',1,...
            'normalize',1,'scale',1,'rho_x',1)
        variable x(n)
        maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
        sum(x) == B
        x >= 0
        if (save_results)
            output = evalc('cvx_end')
        else
            output='';
            cvx_end
        end
        toc
        
        scs_direct.x{i} = x;
        scs_direct.x_viol{i} = min(x);
        scs_direct.budget_viol{i} = abs(1-sum(x));
        scs_direct.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
        scs_direct.output{i} = output;
        
        if (save_results); save('data/portfolio_scs_direct', 'scs_direct'); end
    end
    if run_scs_indirect
        %%
        tic
        cvx_begin
        cvx_solver scs
        cvx_solver_settings('use_indirect',1,'eps',1e-3,'scale',1,'cg_rate',1.5)
        variable x(n)
        maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
        sum(x) == B
        x >= 0
        if (save_results)
            output = evalc('cvx_end')
        else
            output='';
            cvx_end
        end
        toc
        
        scs_indirect.x{i} = x;
        scs_indirect.x_viol{i} = min(x);
        scs_indirect.budget_viol{i} = abs(1-sum(x));
        scs_indirect.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
        scs_indirect.output{i} = output;
        
        
        if (save_results); save('data/portfolio_scs_indirect', 'scs_indirect'); end
    end
    %%
    if run_cvx
        try
            tic
            cvx_begin
            cvx_solver(cvx_use_solver)
            variable x(n)
            maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
            sum(x) == B
            x >= 0
            if (save_results)
                output = evalc('cvx_end')
            else
                output='';
                cvx_end
            end
            toc
            
            cvx.x{i} = x;
            cvx.x_viol{i} = min(x);
            cvx.budget_viol{i} = abs(1-sum(x));
            cvx.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
            cvx.output{i} = output;
            cvx.err{i} = 0;
            
        catch err
            err
            cvx.err{i} = err;
        end
        
        if (save_results); save(sprintf('data/portfolio_cvx_%s',cvx_use_solver), 'cvx'); end
        
    end
end
