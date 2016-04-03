% unit test for linear system solve
n = 100;
m = 200;
num_tests = 50;
for i = 1:num_tests
    % randomly generate data
    data.A = randn(m, n);
    data.b = randn(m, 1);
    data.c = randn(n, 1);
    rho_x = 0.1234; % arbitrary
    % form I + Q
    IQ = sparse([rho_x * eye(n) data.A' data.c;
        -data.A eye(m) data.b;
        -data.c' -data.b' 1]);
    % random right hand side
    u = randn(n+m+1,1);
    v = randn(n+m+1,1);
    % true solution to test against (I + Q)^{-1} * (u+v)
    ut = u+v;
    ut(1:n) = rho_x*ut(1:n);
    true_solution = IQ \ ut;
    clear ut IQ % ensure no leakage
    
    % set up for SCS
    h = [data.c;data.b];
    cg_rate = 2;
    % factorize and form pre-conditioner
    W=sparse([rho_x*speye(n) data.A';data.A -speye(m)]);
    [work.L,work.d,work.P] = ldl(W,'vector');
    work.M = 1 ./ diag(rho_x*speye(n) + data.A'*data.A); % pre-conditioner
    
    % solve both direct and indirect
    for use_indirect = 0:1
        % solve for g, as in paper
        [g, gTh, ~] = solve_for_g(work, data, h, n, m, rho_x, use_indirect ,cg_rate,0);
        % solve full system using cached quantities
        [scs_solution, ~] = project_lin_sys(work, data, n, m, u, v, rho_x, -1, use_indirect, cg_rate, 0, h, g, gTh);
        err = norm(scs_solution - true_solution)
        if err > 1e-8
            fprintf('XXXX incorrect solution! XXXX\n')
        end
    end
end