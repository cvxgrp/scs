function [d,S,U,sig] = ls_broyden(Rx,Rxold,tau,Sk,Yk,S,U,sig,iter,optsDir)
im = mod(iter-1,optsDir.memory);

YSk = Yk'*Sk;
qf = -tau*(Sk'*Rxold);
if YSk<optsDir.delta*abs(qf)
    theta = (1-sgn(qf)*optsDir.delta)*qf/(qf-YSk);
    Yk = theta*Yk-tau*(1-theta)*Rxold;
end

if im == 0 % H = Id
    d = - Rx;
    S = [];U = [];
else
    HYk = Yk;
    for j=1:im-1
        HYk = HYk + U(:,j)*(S(:,j)'*HYk);
    end
    S = [S Sk];
    U = [U (Sk-HYk)/(Sk'*HYk)];
    d = -Rx;
    for j=1:im
        d = d + U(:,j)*(S(:,j)'*d);
    end
end
function x=sgn(x)
if x == 0
    x = 1;
    return
end
x = sign(x);
