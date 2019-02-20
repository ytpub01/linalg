% yasser@ucla.edu math 270c hw 1 CG Steepest Descent Spring 2011
%
function [x,flag,relres,iter,resvec] = mycg_steepest(A,b,tol,maxiter)

k=0;
x=zeros(size(b));
r=b;
ro=norm(r,2)^2;
resvec=zeros(maxiter,1);
p=r;

while (sqrt(ro)/norm(b,2) >= tol && k < maxiter)
    k=k+1;
    Ap=A*p;
    if abs(p'*Ap)/norm(b,2) < tol
        error('Divide by zero');
    end
    alpha=ro/(p'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    if ro/norm(b,2) < tol
        error('Divide by zero');
    end
    beta=norm(r,2)^2/ro;
    ro=norm(r,2)^2;
    resvec(k,1)=norm(r,2);
    p=r+beta*p;
end
if norm(r,2)/norm(b,2) <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/norm(b,2);
iter=k;