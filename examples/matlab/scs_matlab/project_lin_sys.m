function [ut, itn] = project_lin_sys( work, data, n, m, u, v, rho_x, i, use_indirect, cg_rate, extra_verbose, h, g, gTh)
    ut = u+v;
    ut(1:n) = rho_x*ut(1:n);
    ut(1:n+m) = ut(1:n+m) - ut(end)*h;
    ut(1:n+m) = ut(1:n+m) - h*((g'*ut(1:n+m))/(gTh+1));
    warm_start = u(1:n+m);
    ut(n+1:end-1) = -ut(n+1:end-1);
    [ut(1:n+m), itn] = solve_lin_sys(work, data, ut(1:n+m), n, m, warm_start, rho_x, i, use_indirect, cg_rate, extra_verbose);
    ut(end) = (ut(end) + h'*ut(1:n+m));
end