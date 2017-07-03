clear
rng('default');
rng(1);

density = 0.1;

p = 200; % features
q = 100*p; % total samples
w_true = sprandn(p,1,0.2);

X_tmp = 3*sprandn(p, q, density);
ips = -w_true'*X_tmp;
ps = (exp(ips)./(1 + exp(ips)))';
labels = 2*(rand(q,1) < ps) - 1;

X_pos = X_tmp(:,labels==1);
X_neg = X_tmp(:,labels==-1);

X = [X_pos -X_neg]; % include labels with data

lam = 2; 
%%
tic
cvx_begin
    cvx_solver scs
    cvx_solver_settings('eps', 1e-3,...
        'do_super_scs', 1,...
        'direction', 100,...
        'ls', 5,...
        'memory', 50,...
        'verbose', 2)
    variable w(p)
    minimize(sum(log_sum_exp([sparse(1,q);w'*X])) + lam * norm(w,1))        
cvx_end
toc