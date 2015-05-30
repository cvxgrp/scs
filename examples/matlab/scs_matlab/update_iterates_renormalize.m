function [u,v] = update_iterates_renormalize(u, v, D, E, sc_c, sc_b, D_old, E_old, sc_b_old, sc_c_old)
n = size(E,1);
m = size(D,1);
u(1:n) = (u(1:n) ./ (E_old * sc_b_old)) .* (E * sc_b);
u(n+1:n+m) = (u(n+1:n+m) ./ (D_old * sc_c_old)) .* (D * sc_c);
v(n+1:n+m) = (v(n+1:n+m) .* (D_old / (sc_b_old))) ./ (D / (sc_b));
v(end) = (v(end) / (sc_b_old * sc_c_old)) * (sc_b * sc_c);