function z = proj_sdp(z,n)
if n==0
    return;
elseif n==1
    z = max(z,0);
    return;
end

% expand to full size matrix
b = tril(ones(n));
b(b == 1) = z;
z = b;
z = (z + z');
z = z - diag(diag(z)) / 2;

% rescale so projection works, and matrix norm preserved
% see http://www.seas.ucla.edu/~vandenbe/publications/mlbook.pdf pg 3
% scale diags by sqrt(2)
z(eye(n) == 1) = z(eye(n) == 1) .* sqrt(2);

[V,S] = eig(z);
S = diag(S);

idx = find(S>0);
V = V(:,idx);
S = S(idx);
z = V*diag(S)*V';

% scale diags by 1/sqrt(2)
z(eye(n) == 1) = z(eye(n) == 1) ./ sqrt(2);

z = z(tril(ones(n)) == 1);
end