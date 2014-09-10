function [x, i] = pcg_scs(A,b,x,M,rho_x,max_iters,tol,verbose)
% M is inverse preconditioner
r=b-(rho_x*x+A'*(A*x));
z = M.*r;
p = z;
ip = r'*z;
for i=1:max_iters
    Ap=rho_x*p + A'*(A*p);
    alpha= ip/(p'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    if norm(r)<tol
        if verbose
            fprintf('CG took %i iterations to converge, resisdual %4f <= tolerance %4f\n',i,resid,tol)
        end
        return;
    end
    z = M.*r;
    ipold = ip;
    ip = z'*r;
    beta =  ip / ipold;
    p=z+beta*p;
end
if verbose
    fprintf('CG did not converge within %i iterations, resisdual %4f > tolerance %4f\n',max_iters,resid,tol)
end
end
