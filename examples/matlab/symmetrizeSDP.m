function z = symmetrizeSDP(z,K)
l = K.f + K.l;
for i=1:length(K.q)
    l = l + K.q(i);
end
for i=1:length(K.s)
    V = reshape(z(l+1:l+K.s(i)^2),K.s(i),K.s(i));
    V = (V+V')/2;
    z(l+1:l+K.s(i)^2) = V(:);
    l=l+K.s(i)^2;
end
