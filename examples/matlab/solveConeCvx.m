function [x,y,s]=solveConeCvx(data,K,cs)
if (K.ep>0 || K.ed>0)
    x = nan;
    y = nan;
    s = nan;
    disp('cvx cannot solve ECPs');
    return;
end
n = length(data.c);
m = length(data.b);
% can NOT solve ECPs
cvx_begin
cvx_solver(cs)
variables x(n) s(m)
dual variable y
minimize(data.c'*x)
y:data.A*x + s == data.b
s(1:K.f)==0
l = K.f;
s(l+1:l+K.l) >= 0
l = l + K.l;
for i=1:length(K.q)
    s(l+1) >= norm(s(l+2:l + K.q(i)));
    l = l + K.q(i);
end
for i=1:length(K.s)
    reshape(s(l+1:l + (K.s(i))^2),K.s(i),K.s(i)) == semidefinite(K.s(i));
    l = l + (K.s(i))^2;
end
cvx_end
end
