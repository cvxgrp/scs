function [data, w] = unnormalize_data(data, scale, w)
data.A = sparse(diag(1./w.D))*data.A*sparse(diag(1./w.E) / scale);

data.b = data.b .* w.D / (w.sc_b * scale);
data.c = data.c .* w.E / (w.sc_c * scale);
end