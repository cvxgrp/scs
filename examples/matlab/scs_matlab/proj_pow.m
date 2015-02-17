function [x, iter] = proj_pow(in, a)

CONE_MAX_ITERS = 20;
CONE_TOL = 1e-8;

iter = 0;
x = in(1); y = in(2); z = in(3);
% v in cl(K_a)
if(x >=0 && y>=0 && (x^a) * (y^(1-a)) >= abs(z));
    x = in;
    return
end

% -v in K_a^*
u = -in(1); v = -in(2); w = -in(3);
if (u>=0 && v>=0 && (u^a) * (v^(1-a)) >= abs(w) * (a^a) * ((1-a)^(1-a)));
    x = zeros(3,1);
    return
end

xh = in(1);
yh = in(2);
zh = in(3);
rh = abs(zh);
r = rh / 2;
for iter=1:CONE_MAX_ITERS;
    x = calcX(r, xh, rh, a);
    y = calcX(r, yh, rh, 1-a);
    
    f = calcF(x,y,r,a);
    if abs(f) < CONE_TOL
        break
    end
    
    dxdr = calcdxdr(x,xh,rh,r,a);
    dydr = calcdxdr(y,yh,rh,r,(1-a));
    fp = calcFp(x,y,dxdr,dydr,a);
    
    r = min(max(r - f/fp,0), rh);
end
z = sign(zh) * r;
x = [x;y;z];
end

function x = calcX(r, xh, rh, a)
x = max(0.5 * (xh + sqrt(xh*xh + 4 * a * (rh - r) * r)), 1e-12);
end

function dx = calcdxdr(x,xh,rh,r, a)
dx = a * (rh - 2*r) / (2*x - xh);
end

function f = calcF(x,y,r,a)
f = (x^a) * (y^(1-a)) - r;
end

function fp = calcFp(x,y,dxdr,dydr,a)
fp = (x^a) * (y^(1-a)) * (a * dxdr / x + (1-a) * dydr / y) - 1;
end