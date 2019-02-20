% yasser@ucla.edu math 270c hw 1 CG Lanczos Spring 2011
%
function [x,flag,relres,k,resvec] = mycg_lanczos(A,b,tol,maxiter)

%compute q1
q=b/norm(b,2);
x=zeros(size(b));
residual=b-A*x;
resvec=zeros(maxiter,1);
Aq=A*q;                            %compute h1,1
alpha=q'*Aq;
d=alpha;
if (abs(d)/norm(b,2) > tol^2)
    p=norm(b,2)/d;
else
    error('Divide by zero');
end
c=q;
x=c*p;
v=Aq-alpha*q; %v1
beta=norm(v,2);

k=1;
while (k < maxiter && abs(beta) > tol)
    k=k+1;
    resvec(k,1)=norm(residual,2);
    qprev=q;
    q=v/norm(v,2);
    Aq=A*q;
    alpha=q'*Aq;
    if (abs(d)/norm(b,2) > tol^2)
        miu=beta/d;
    else
        error('Divide by zero');
    end
    d=alpha-miu*beta;
    if (abs(d)/norm(b,2) > tol^2)
        p=-p*beta/d;
    else
        error('Divide by zero');
    end
    c=q-miu*c;
    x=x+c*p;
    v=Aq-alpha*q-beta*qprev;
    beta=norm(v,2);
    
    residual=b-A*x;
    if norm(residual-resvec(k,1),2)/norm(b,2) < tol
        error('CG stagnated, exiting');
    end
end
if norm(residual,2)/norm(b,2) <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/norm(b,2);