% yasser@ucla.edu math 270c hw 1 Gauss-Siedel Spring 2011
%
function [x,flag,relres,k,resvec] = mysiedel(A,b,tol,maxiter)

x=zeros(size(b));
resvec=zeros(maxiter,1);
r=b;

k=0;
while (k < maxiter && norm(r,2) >= tol)
    k=k+1;
    resvec(k,1)=norm(r,2);
    N=size(A,1);
    for i=1:N
        x(i)=x(i)+r(i)/A(i,i);
        r=b-A*x;
    end
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
