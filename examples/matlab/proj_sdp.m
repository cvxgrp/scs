function z = proj_sdp(z,n)
if isempty(z)
    return;
elseif length(z)==1
    z = max(z,0);
    return;
elseif length(z)==2
    z = proj_sdp2(z);
    return;
end
z = reshape(z,n,n);
zs=(z+z')/2;

[V,S] = eig(zs);
S = diag(S);

idx = find(S>0);
V = V(:,idx);
S = S(idx);
z = V*diag(S)*V';

z = z(:);
end