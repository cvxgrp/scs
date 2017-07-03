clear;
rng('default');
rng(1);

d = 200;
p = 5;
A = sprandn(p,d,0.3,0.0001);
S = full(A'*A);

lam = 2;

cvx_begin sdp
    cvx_solver scs
    cvx_solver_settings('eps', 1e-3,...
        'verbose', 1,...
        'do_super_scs', 0, ...
        'direction', 100, ...
        'memory', 100);
    variable X(d,d) symmetric
    minimize(-trace(S*X) + lam*norm(X,1))
    trace(X)==1
    X>=0
cvx_end  
