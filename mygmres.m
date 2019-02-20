% yasser@ucla.edu math 270c hw 1 GMRES Spring 2011
%
function [x,flag,relres,iter,resvec] = mygmres(A,b,tol,maxiter)

%compute q1
q=b/norm(b,2);
Q=q;
if maxiter > length(b)
    warning(['MAXIT is ' num2str(maxiter) .... 
        ' but it should be bounded by SIZE(A,1).' ...
         ' Setting MAXIT to ' num2str(length(b)) '.'])
    maxiter=min(maxiter,length(b));
end
Aq=A*q;                            %compute h1,1
H=q'*Aq;
R=[];
v=Aq-Q*H(:,1); %v1
bhat=norm(b,2);
c_givens=bhat;
x=zeros(size(b));
residual=b-A*x;
resvec=zeros(maxiter,1);

k=0;
while (k < maxiter && norm(v,2)/norm(b,2) >= tol...
        && (norm(residual,2)/norm(b,2) >= tol))
    k=k+1;
    resvec(k,1)=norm(residual,2);
    q=v/norm(v,2); %new search direction
    %update H
    H=[H;zeros(1,k)];
    H(k+1,k)=norm(v,2);
    bhat=[norm(b,2);zeros(k,1)];
    %don't start over with c_givens, get from last
    c_givens=[c_givens;0];
%     r=zeros(k+1,1);
    r=H(:,k);
    
    % least squares
    % Q R of H
    for i=1:k-1
        nrm=norm([c(i,1);s(i,1)],2);
        G=nrm*eye(k+1,k+1);
        G(i,i)=c(i,1);
        G(i+1,i+1)=c(i,1);
        G(i,i+1)=s(i,1);
        G(i+1,i)=-s(i,1);
        G=sparse(G);
        %get last column of R
        if nrm/norm(b,2) < tol
            error('Divide by zero');
        end;
        r=G'*r/nrm; %forces r(i+1)=0
    end;
    
    %last column of R
    c(k,1)=r(k,1);
    s(k,1)=-r(k+1,1);
    nrm=norm([c(k,1);s(k,1)]);
    G=nrm*eye(k+1,k+1);
    G(k,k)=c(k,1);
    G(k+1,k+1)=c(k,1);
    G(k,k+1)=s(k,1);
    G(k+1,k)=-s(k,1);
    G=sparse(G);
    
    %get last column of R
    if nrm/norm(b,2) < tol
        error('Divide by zero');
    end;
    r=G'*r/nrm; %forces r(k+1)=0
    c_givens=G'*c_givens/nrm;
    
    %set new R
    if k>1
        R=[R; zeros(1,k-1)];
    end;
    R=[R r(1:end-1,1)];
    
    %get lambda
    lambda=zeros(k,1);
    for i=0:k-1
        lambda(k-i,1)=c_givens(k-i,1);
        for j=k-i+1:k
            lambda(k-i,1)=lambda(k-i,1)-R(k-i,j)*lambda(j,1);
        end
        lambda(k-i,1)=lambda(k-i,1)/R(k-i,k-i);
    end
    %set x
    x=Q*lambda;
    %end least squares
    
    %set Qk+1 with qk+1
    Q=[Q q];
    %set Hk+1 except last row
    Aq=A*q;
    H=[H Q'*Aq];
    v=Aq-Q*H(:,k+1); %vk+1
    residual=b-A*x;
    if norm(residual-resvec(k,1),2)/norm(b,2) < tol
        error('GMRES stagnated, exiting');
    end
end
if norm(residual,2)/norm(b,2) <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/norm(b,2);
iter=k;