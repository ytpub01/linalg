% yasser@ucla.edu math 270c hw 1 GMRES Spring 2011
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
    xk=zeros(N-1,1);
    
    e=-1*ones(1,N-2);
    d=2*ones(1,N-1);
    A3=gallery('tridiag',e,d,e);
    A3(N-1,N-2)=-dx;
    A3(N-1,N-1)=dx;
    A3=(1/dx)^2*A3;
    b3=circshift(f,-1);
    b3(end-1:end)=[];
    b3(1)=b3(1)+a/(dx^2);
    b3(end)=b;
    
    [u3 flag relres iter rkperN]=mygmres(A3,b3,tol,maxiter);
    
    u3=[a;u3;a];
    eN=[eN norm(uexact-[a;A3\b3;a],Inf)];
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
    plot(x,uexact,'-b',x,[0;A3\b3;0],'-.g',x,u3,'-.r');
    xlabel('x');
    ylabel('y');
    title(['Solution to the PDE at iteration k=' num2str(iter(1))...
        ' for N=' num2str(N)]);
    legend('exact','discretized','GMRES','Location','SouthWest');
end

slope=polyfit(log10(Nvalues),log10(eN),1);
figure;
plot(log10(Nvalues),log10(eN),'-.or');
xlabel('log10 ( N )');
ylabel('log10 ( max ( e^N ) )');
title(['Discretization error of the PDE for A3, slope = ' num2str(slope(1))]);
figure;
plot(1:rkiter(1),log10(rk(1,1:rkiter(1))),1:rkiter(2),...
    log10(rk(2,1:rkiter(2))),1:rkiter(3),log10(rk(3,1:rkiter(3))),...
    1:rkiter(4),log10(rk(4,1:rkiter(4))));
xlabel('k');
ylabel('log10 ( norm ( r_k ) )');
title('L2 norm of the GMRES residual at each iteration');
legend(num2str(Nvalues(1)),num2str(Nvalues(2)),num2str(Nvalues(3)),...
    num2str(Nvalues(4)),'Location','SouthWest');
