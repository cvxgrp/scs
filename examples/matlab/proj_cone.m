function z = proj_cone(z,c)
free_len = c.f;
lp_len = c.l;
k_soc = length(c.q);
q = c.q;
s = c.s;
ssize = length(c.s);
% free/zero cone
z(1:free_len) = 0;
% lp cone
z(free_len+1:lp_len+free_len) = max(z(free_len+1:lp_len+free_len),0);
% SOCs
idx=lp_len+free_len;
for i=1:k_soc
    z(idx+1:idx+q(i)) = proj_soc(z(idx+1:idx+q(i)));
    idx=idx+q(i);
end
%SDCs
for i=1:ssize
    z(idx+1:idx+s(i)^2) = proj_sdp(z(idx+1:idx+s(i)^2),s(i));
    idx=idx+s(i)^2;
end
%Exp primal
for i=1:c.ep
    z(idx+1:idx+3) = project_exp_bisection(z(idx+1:idx+3));
    idx = idx+3;
end
%Exp dual
for i=1:c.ed
    z(idx+1:idx+3) = z(idx+1:idx+3) + project_exp_bisection(-z(idx+1:idx+3));
    idx = idx+3;
end

end

function z = proj_soc(tt)
if isempty(tt)
    z=[];
    return;
elseif length(tt)==1
    z = max(tt,0);
    return;
end
v1=tt(1);v2=tt(2:end);
if norm(v2)<=-v1
    v2=zeros(length(v2),1);
    v1=0;
elseif norm(v2) > abs(v1)
    v2=0.5*(1+v1/norm(v2))*v2;
    v1=norm(v2);
end
z=[v1;v2];
end

function z = proj_sdp(z,n)
if isempty(z)
    z=[];
    return;
elseif length(z)==1
    z = max(z,0);
    return;
end

z = reshape(z,n,n);
zs=(z+z')/2;

%ii
[V,S] = eig(zs);
S = diag(S);
num_pos = sum(S>0);
num_neg = sum(S<0);
if (num_pos < num_neg)
    positive = true;
    idx = find(S>0);
    V = V(:,idx);
    S = S(idx);
else
    positive = false;
    idx = find(S<0);
    V = V(:,idx);
    S = S(idx);
end

if (positive)
    T = S;
    T(T<0) = 0;
    z = V*diag(T)*V';
else
    T = S;
    T(T>0) = 0;
    z = zs - V*diag(T)*V';
end
z = z(:);
end

function x = project_exp_bisection(v)
r = v(1); s = v(2); t = v(3);
% v in cl(Kexp)
if( (s.*exp(r./s) <= t && s > 0) || (r <= 0 && s == 0 && t >= 0) );
    x = v;
    return
end
x = zeros(3,1);
% -v in Kexp^*
if ( (-r < 0 && r.*exp(s./r) <= -exp(1).*t) || (-r == 0 && -s >= 0 && -t >= 0) );
    return
end

% special case with analytical solution
if(r < 0 && s < 0);
    x = v;
    x(2) = 0;
    x(3) = max(v(3),0);
    return
end

x = v;
[ub,lb] = getRhoUb(v);
for iter=1:1e2;
    rho = (ub + lb)/2;
    [g,x] = calcGrad(v,rho,x);
    if (g > 0)
        lb = rho;
    else
        ub = rho;
    end
    if (ub - lb < 1e-6)
        break
    end
end
end

function [ub,lb] = getRhoUb(v)
lb = 0;
rho = 2^(-3);
[g,z] = calcGrad(v,rho,v);
while (g>0)
    lb = rho;
    rho = rho*2;
    [g,z] = calcGrad(v,rho,z);
end
ub = rho;
end

function [g,x] = calcGrad(v,rho,warm_start)
x = solve_with_rho(v,rho,warm_start(3));
if (x(2)==0)
    g = x(1);
else
    g = (x(1) + x(2)*log(x(2)/x(3)));
end
end


function x = solve_with_rho(v,rho,w)
x = zeros(3,1);
x(3) = newton_exp_onz(rho,v(2),v(3),w);
x(2) = (1/rho)*(x(3) - v(3))*x(3);
x(1) = v(1) - rho;
end


function z = newton_exp_onz(rho, y_hat, z_hat,w)
t = max(max(w - z_hat, -z_hat),1e-6);
for iter=1:1e2;
    f = (1/rho^2)*t*(t + z_hat) - y_hat/rho + log(t/rho) + 1;
    fp = (1/rho^2)*(2*t + z_hat) + 1/t;
    
    t = t - f/fp;
    if (t <= -z_hat)
        t = -z_hat;
        break;
    elseif (t <= 0)
        t = 0;
        break;
    elseif (abs(f)<1e-6)
        break;
    end
end
z = t + z_hat;
end