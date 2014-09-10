function [y, itn] = solve_lin_sys(w, data, rhs, n, m, warm_start, rho_x, iter, use_indirect, cg_rate, extra_verbose)
% assumes matrix [rho_x*I A';A -I]
BEST_TOL = 1e-9;
if use_indirect
    if (iter < 0)
        tol = BEST_TOL;
    else
        tol = 0.1 * (iter+1)^(-cg_rate);
    end
    tol = max( BEST_TOL, norm(rhs(1:n)) * tol);
    y = zeros(n+m,1);
    [y(1:n), itn] = pcg_scs(data.A,rhs(1:n)+data.A'*rhs(n+1:n+m),warm_start(1:n),w.M,rho_x,n,tol,extra_verbose);
    y(n+1:n+m) = -rhs(n+1:n+m) + data.A*y(1:n);
else
    % SOLVE DIRECTLY:
    y = (w.L'\(w.d\(w.L\rhs(w.P))));y(w.P)=y;
    itn = -1;
end
end