% yasser@ucla.edu math 270c hw 2 MINRES2 Spring 2011
%
function [x,flag,relres,iter,resvec] = myminres2(A,b,tol,maxiter)

resvec=zeros(maxiter,1);
normb=norm(b,2);
x=zeros(size(b));                                       %x0=0

k=1;
q=b/normb;                                              %q1
Aq=A*q;                                                 %Aq1
alpha=q'*Aq;                                            %alpha1
gamma_bar=alpha;                                        %gamma_bar1:alpha1
v=Aq-alpha*q;                                           %v2

%for minres2
beta=norm(v,2);                                         %beta2
c=gamma_bar;                                            %gamma_bar1 so c1
s=-beta;                                                %beta2 so s1
gamma=norm([c s],2);                                    %gamma1
tau=normb;                                              %tau1
if abs(gamma)/normb < tol
    error('Divide by zero');
end
t=tau*c/gamma;                                          %t1:c1
m=1/gamma*q;                                            %m1:gamma1,q1
x=x+t*m;                                                %x1

residual=b-A*x;                                         %r1
resvec(k,1)=norm(residual,2);

m_prev=m;                                               %m1

k=2;
qprev=q;                                                %q1
if beta/normb < tol
    error('Divide by zero');                            %beta2
end
q=v/beta;                                               %q2:v2,beta2
Aq=A*q;                                                 %Aq2
alpha=q'*Aq;                                            %alpha2
gamma_bar=(s*beta+c*alpha)/gamma;                       %gamma_bar2
delta=(c*beta-s*alpha)/gamma;                           %delta2
cprev=c/gamma;                                          %c1
sprev=s/gamma;                                          %s1
v=Aq-alpha*q-beta*qprev;                         %v3:Aq2,alpha2,q2,beta2,q1

%for minres2
beta=norm(v,2);                                         %beta3
c=gamma_bar;                                            %gamma_bar2 so c2
s=-beta;                                                %beta3 so s2
gamma=norm([c s],2);                                    %gamma2
tau=sprev*tau;
if abs(gamma)/normb < tol
    error('Divide by zero');
end
t=tau*c/gamma;                                          %t2:s1,c2,tau1
m=1/gamma*(q-delta*m_prev);                         %m2:gamma2,q2,delta2,m1
x=x+t*m;                                                %x2

residual=b-A*x;                                         %r2
if norm(residual-resvec(k-1,1),2)/normb < tol
    error('MINRES1 stagnated, exiting');
end
resvec(k,1)=norm(residual,2);

m_prevprev=m_prev;                                       %m1
m_prev=m;                                               %m2

while (k < maxiter && norm(v,2)/normb >= tol...
        && (norm(residual,2)/normb >= tol))
    k=k+1;                                              %k=3
    qprev=q;                                            %q2
    if beta/normb < tol
        error('Divide by zero');                        %beta3
    end
    q=v/beta;                                           %q3
    Aq=A*q;                                             %Aq3
    alpha=q'*Aq;                                        %alpha3
    gamma_bar=(s*cprev*beta+c*alpha)/gamma;%gamma_br3:s2,c1,beta3,c2,alpha3
    delta=(c*cprev*beta-s*alpha)/gamma;      %delta3:c2,c1,beta3,s2,alpha3
    epsilon=-sprev*beta;                                %epsilon3:s1,beta3
    cprev=c/gamma;                                      %c2
    sprev=s/gamma;                                      %s2
    v=Aq-alpha*q-beta*qprev;                    %v4:Aq3,alpha3,q3,beta3,q2
    
    %new for minres2
    beta=norm(v,2);                                     %beta4
    c=gamma_bar;                                        %gamma_bar3 so c3
    s=-beta;                                            %beta4 so s3
    gamma=norm([c s],2);                                %gamma3
    tau=sprev*tau;                                      %tau3
    if abs(gamma)/normb < tol
        error('Divide by zero');
    end
    t=tau*c/gamma;                                      %t3:s2,c3,tau3
    m=1/gamma*(q-delta*m_prev-epsilon*m_prevprev);%m3:gmma3,q3,delta3,m2,m1
    x=x+t*m;                                            %x3

    residual=b-A*x; 
    if norm(residual-resvec(k-1,1),2)/normb < tol
        error('MINRES1 stagnated, exiting');
    end
    resvec(k,1)=norm(residual,2);
    
    m_prevprev=m_prev;                                  %m2
    m_prev=m;                                           %m3

end
if norm(residual,2)/normb <= tol
    flag=0;
else
    flag=1;
end
resvec=resvec(1:k);
relres=resvec/normb;
iter=k;