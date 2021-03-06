% yasser@ucla.edu math 270c hw 2 MINRES2 Spring 2011
%

close all
clear all
format longG
a=0;
b=2*pi;
tol=1e-10;
maxiter=64;
Nvalues=[65,129,257,513];
eN=[];
rk=[];
rkiter=[];

for N=Nvalues
    dx=1/N;
    x=0:dx:1;
    uexact=sin(2*pi*x)';
    f=4*pi^2*sin(2*pi*x)';
    f_bar=(1/dx)*4*pi*sin(pi*dx)*sin(2*pi*x+pi*dx)';
    
    offdiag=-1*ones(1,N);
    diag=2*ones(1,N+1);
    A2=gallery('tridiag',offdiag,diag,offdiag);
    A2(1,1)=1;
    A2(N+1,1)=dx;
    A2(1,N+1)=dx;
    A2(N,N)=1;
    A2(N+1,N)=0;
    A2(N,N+1)=0;
    A2(N+1,N+1)=0;
    A2=1/dx*A2;
    b2=(circshift(f_bar,1)+f_bar)/2*dx;
    b2(1,1)=f_bar(1,1)/2*dx;
    b2(end-1)=f_bar(N-1,1)/2*dx+b;
    b2(end)=a;
    
    [u3 flag relres iter rkperN]=myminres2(A2,b2,tol,maxiter);
    
    x=x(1:end-1);   %last element not part of the solution
    uexact=uexact(1:end-1);
    u3=u3(1:end-1);
    udiscret=A2\b2;
    udiscret=udiscret(1:end-1);
    eN=[eN norm(uexact-udiscret,Inf)];
    if maxiter > size(rkperN,1)
        fill=zeros(1,maxiter-size(rkperN,1));
        fill(1,:)=rkperN(end);
        rk=[rk;rkperN' fill];
    elseif maxiter < size(rkperN,1)
        rk=[rk;rkperN(1:maxiter)'];
    else
        rk=[rk;rkperN'];
    end
    rkiter=[rkiter;iter];
    figure;
    plot(x,uexact,'-b',x,udiscret,'-.g',x,u3,'-.r');
    xlabel('x');
    ylabel('y');
    title(['Solution to the PDE at iteration k=' num2str(iter)...
        ' for N=' num2str(N)]);
    legend('exact','discretized','MINRES2','Location','SouthWest');
end

slope=polyfit(log10(Nvalues),log10(eN),1);
figure;
plot(log10(Nvalues),log10(eN),'-.or');
xlabel('log10 ( N )');
ylabel('log10 ( max ( e^N ) )');
title(['Discretization error of the PDE for A2, slope = ' num2str(slope(1))]);
figure;
plot(1:rkiter(1),log10(rk(1,1:rkiter(1))),1:rkiter(2),...
    log10(rk(2,1:rkiter(2))),1:rkiter(3),log10(rk(3,1:rkiter(3))),...
    1:rkiter(4),log10(rk(4,1:rkiter(4))));
xlabel('k');
ylabel('log10 ( norm ( r_k ) )');
title('L2 norm of the MINRES2 residual at each iteration');
legend(num2str(Nvalues(1)),num2str(Nvalues(2)),num2str(Nvalues(3)),...
    num2str(Nvalues(4)),'Location','SouthWest');
