% yasser@ucla.edu math 270c hw 2 GMRES with restarts Spring 2011
%
function [x,flag,relres,iter,resvec] = mygmres_restart(A,b,m,tol,maxiter)

k=0;
rhs=b;
x=zeros(size(b));
residual=[];
div=floor(maxiter/m);
for i=1:div
    [y,flag,relres,iter,resvec] = mygmres(A,rhs,tol,m);
    residual=[residual;resvec];
    k=k+iter;
    x=x+y;
    rhs=b-A*x;
end
remainder=mod(maxiter,m);
if remainder > 0
    [y,flag,relres,iter,resvec] = mygmres(A,rhs,tol,remainder);
    residual=[residual;resvec];
    k=k+iter;
    x=x+y;
end
resvec=residual;
relres=residual/norm(b,2);
iter=[k iter];
