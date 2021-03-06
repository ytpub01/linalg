% yasser@ucla.edu math 270c hw 2 BICG Steepest Descent Spring 2011
%
function [x,flag,relres1,relres2,iter,resvec1,resvec2] = mybicg_steepest(A,b,tol,maxiter)

k=0;
x=zeros(size(b));
r=b;
s=b;
residual1=norm(r,2)^2;
resvec1=zeros(maxiter,1);
residual2=norm(s,2)^2;
resvec2=zeros(maxiter,1);
p=r;
q=r;
ro=s'*r;

while (sqrt(residual1) >= tol && sqrt(residual2) >= tol && k < maxiter)
    k=k+1;
    Ap=A*p;
    if abs(q'*Ap)/norm(b,2) < tol
        error('Divide by zero');
    end
    alpha=ro/(q'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    residual1=norm(r,2);
    resvec1(k,1)=residual1;
    Aq=A*q;
    s=s-alpha*Aq;
    residual2=norm(s,2);
    resvec2(k,1)=residual2;
    if ro/norm(b,2) < tol
        error('Divide by zero');
    end
    beta=s'*r/ro;
    ro=s'*r;
    p=r+beta*p;
    q=s+beta*q;
end
if (residual1/norm(b,2) <= tol && residual2/norm(b,2) <= tol)
    flag=0;
else
    flag=1;
end
resvec1=resvec1(1:k);
resvec2=resvec2(1:k);
relres1=resvec1/norm(b,2);
relres2=resvec2/norm(b,2);
iter=k;
