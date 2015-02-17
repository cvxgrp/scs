function K = validateCone(K)
if (~isfield(K, 'f'))
    K.f = 0;
end
if (~isfield(K, 'l'))
    K.l = 0;
end
if (~isfield(K, 'ep'))
    K.ep = 0;
end
if (~isfield(K, 'ed'))
    K.ed = 0;
end
if (~isfield(K, 'q'))
    K.q = [];
end
if (~isfield(K, 's'))
    K.s = [];
end
if (~isfield(K, 'p'))
    K.p = [];
end
end
