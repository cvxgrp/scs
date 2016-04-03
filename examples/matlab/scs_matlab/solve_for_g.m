function [g, gTh, itn] = solve_for_g(work,data,h,n,m,rho_x,use_indirect,cg_rate,extra_verbose)
    [g, itn] = solve_lin_sys(work,data,h,n,m,zeros(n,1),rho_x,-1,use_indirect,cg_rate,extra_verbose);
    g(n+1:end) = -g(n+1:end);
    gTh = g'*h;
end