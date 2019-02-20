% yasser@ucla.edu math 270c hw 3 PCG Spring 2011
%
function [x,flag,relres,iter,resvec] = mypcg(A,b,tol,maxiter)

k=0;
x=zeros(size(b));
r=b;
r_prev=r;
residual=norm(r,2);
resvec=zeros(maxiter,1);

n=size(A,1);
y=zeros(n,1);
z=zeros(n,1);
p=zeros(n,1);

R=mycholu(A);
L=R';

while (residual >= tol && k < maxiter)
    %solve lower
    residual
    y(1)=r(1)/L(1,1);
    for j=2:n
        y(j)=(r(j)-L(j,1:j-1)*y(1:j-1))/L(j,j);
    end;
    z_prev=z;
    %solve upper
    z(n)=y(n)/R(n,n);
    for j=n-1:-1:1
        z(j)=(y(j)-R(j,j+1:n)*z(j+1:n))/R(j,j);
    end;
    k=k+1;
    if (k==1)
        p=z;
    else
        beta=(r'*z)/(r_prev'*z_prev);
        p=z+beta*p;
        Ap=A*p;
        alpha=(r'*z)/(p'*Ap);
        x=x+alpha*p;
        r_prev=r;
        r=r-alpha*Ap;
        residual=norm(r,2);
        resvec(k,1)=residual;
    end
end
if residual/norm(b,2) <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/norm(b,2);
iter=k;
