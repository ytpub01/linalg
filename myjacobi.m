% yasser@ucla.edu math 270c hw 1 Jacobi Spring 2011
%
function [x,flag,relres,k,resvec] = myjacobi(A,b,tol,maxiter)

x=zeros(size(b));
resvec=zeros(maxiter,1);
r=b;

k=0;
while (k < maxiter && norm(r,2) >= tol)
    k=k+1;
    r=b-A*x;
    resvec(k,1)=norm(r,2);
    N=size(A,1);
    x=x+r./(A(1:N+1:end)');    
    if norm(r-resvec(k,1),2)/norm(b,2) < tol
        error('Jacobi stagnated, exiting');
    end
end
if norm(r,2)/norm(b,2) <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/norm(b,2);
