function x = proj_exp(v)

CONE_MAX_ITERS = 100;
CONE_TOL = 1e-8;

r = v(1); s = v(2); t = v(3);
% v in cl(Kexp)
if( (s*exp(r./s) <= t && s > 0) || (r <= 0 && s == 0 && t >= 0) );
    x = v;
    return
end

% -v in Kexp^*
if ( (-r < 0 && r*exp(s./r) <= -exp(1)*t) || (-r == 0 && -s >= 0 && -t >= 0) );
    x = zeros(3,1);
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
for iter=1:CONE_MAX_ITERS;
    rho = (ub + lb)/2;
    [g,x] = calcGrad(v,rho,x);
    if (g > 0)
        lb = rho;
    else
        ub = rho;
    end
    if (ub - lb < CONE_TOL)
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
t = max(max(w - z_hat, -z_hat),CONE_TOL);
for iter=1:CONE_MAX_ITERS;
    f = (1/rho^2)*t*(t + z_hat) - y_hat/rho + log(t/rho) + 1;
    fp = (1/rho^2)*(2*t + z_hat) + 1/t;
    
    t = t - f/fp;
    if (t <= -z_hat)
        t = -z_hat;
        break;
    elseif (t <= 0)
        t = 0;
        break;
    elseif (abs(f)<CONE_TOL)
        break;
    end
end
z = t + z_hat;
end