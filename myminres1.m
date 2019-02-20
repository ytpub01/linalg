% yasser@ucla.edu math 270c hw 2 MINRES1 Spring 2011
%
function [x,flag,relres,iter,resvec] = myminres1(A,b,tol,maxiter)

x=zeros(size(b));
y=zeros(size(b));
resvec=zeros(maxiter,1);

y=zeros(size(b));
resvec=zeros(maxiter,1);

k=1;            %compute q1
q=b/norm(b,2);
Aq=A*q;
alpha=q'*Aq;
w_bar=q;
if abs(alpha)/norm(b,2) < tol
    error('Divide by zero');
end
z_bar=norm(b,2)/alpha;
x=y+z_bar*w_bar;
gamma_bar=alpha;
residual=b-A*x;
resvec(k,1)=norm(residual,2);
v=Aq-alpha*q;          %v2

k=2;
qprev=q;
if norm(v,2)/norm(b,2) < tol
    error('Divide by zero');
end
q=v/norm(v,2);
Aq=A*q;
alpha=q'*Aq;
beta=norm(v,2);
c=gamma_bar;
s=-beta;
if norm([s c],2)/norm(b,2) < tol
    error('Divide by zero');
end
gamma_bar=(s*beta+c*alpha)/norm([s c],2);
delta=(c*beta-s*alpha)/norm([s c],2);
gamma=norm([c s],2);
if abs(gamma)/norm(b,2) < tol
    error('Divide by zero');
end
z=norm(b,2)/gamma;
if abs(gamma_bar)/norm(b,2) < tol
    error('Divide by zero');
end
z_bar=(-delta*z)/gamma_bar;
w=(c*w_bar-s*q)/norm([c s],2);
w_bar=(s*w_bar+c*q)/norm([c s],2);
y=y+z*w;
x=y+z_bar*w_bar;
z_prev2=z;
z_bar_prev=z_bar;
gamma_bar_prev=gamma_bar;
sprev=s/norm([s c],2);
cprev=c/norm([s c],2);
residual=b-A*x;
if norm(residual-resvec(k-1,1),2)/norm(b,2) < tol
    error('MINRES1 stagnated, exiting');
end
resvec(k,1)=norm(residual,2);
v=Aq-alpha*q-beta*qprev;          %v3

while (k < maxiter && norm(v,2)/norm(b,2) >= tol...
        && (norm(residual,2)/norm(b,2) >= tol))
    k=k+1;
    qprev=q;
    if norm(v,2)/norm(b,2) < tol
        error('Divide by zero');
    end
    q=v/norm(v,2); %new search direction
    Aq=A*q;
    alpha=q'*Aq;
    beta=norm(v,2);
    c=gamma_bar;
    s=-beta;
    if norm([s c],2)/norm(b,2) < tol
        error('Divide by zero');
    end
    gamma_bar=(s*cprev*beta+c*alpha)/norm([s c],2);
    delta=(c*cprev*beta-s*alpha)/norm([s c],2);
    epsilon=-sprev*beta;
    gamma=norm([c s],2);
    if abs(gamma)/norm(b,2) < tol
        error('Divide by zero');
    end
    z=gamma_bar_prev/gamma*z_bar_prev;
    if abs(gamma_bar)/norm(b,2) < tol
        error('Divide by zero');
    end
    z_bar=(-epsilon*z_prev2-delta*z)/gamma_bar;
    w=(c*w_bar-s*q)/norm([c s],2);
    w_bar=(s*w_bar+c*q)/norm([c s],2);
    y=y+z*w;
    x=y+z_bar*w_bar;
    z_prev2=z;
    z_bar_prev=z_bar;
    gamma_bar_prev=gamma_bar;
    sprev=s/norm([s c],2);
    cprev=c/norm([s c],2);
    residual=b-A*x;
    resvec(k,1)=norm(residual,2);
    v=Aq-alpha*q-beta*qprev;          %vk+1
    
    %leave this, ok
    if norm(residual-resvec(k-1,1),2)/norm(b,2) < tol
        error('MINRES1 stagnated, exiting');
    end
    resvec(k,1)=norm(residual,2);

end

% leave this, ok
if norm(residual,2)/norm(b,2) <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/norm(b,2);
iter=k;