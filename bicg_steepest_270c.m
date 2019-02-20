% yasser@ucla.edu math 270c hw 2 BICG Steepest Descent Spring 2011
%

close all
clear all
format longG
a=0;
b=2*pi;
tol=1e-10;
maxiter=33;
Nvalues=[65,129,257,513];
eN=[];
rk1=[];
rk2=[];
rkiter1=[];
rkiter2=[];

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
    
    [u1 flag relres1 relres2 iter rkperN1 rkperN2]=mybicg_steepest(A3,b3,tol,maxiter);

    u1=[a;u1;a];
    eN=[eN norm(uexact-[a;A3\b3;a],Inf)];

    if maxiter > size(rkperN1,1)
        fill=zeros(1,maxiter-size(rkperN1,1));
        fill(1,:)=rkperN1(end);
        rk1=[rk1;rkperN1' fill];
    elseif maxiter < size(rkperN1,1)
        rk1=[rk1;rkperN1(1:maxiter)'];
    else
        rk1=[rk1;rkperN1'];
    end
    rkiter1=[rkiter1;iter];

    if maxiter > size(rkperN2,1)
        fill=zeros(1,maxiter-size(rkperN2,1));
        fill(1,:)=rkperN2(end);
        rk2=[rk2;rkperN2' fill];
    elseif maxiter < size(rkperN2,1)
        rk2=[rk2;rkperN2(1:maxiter)'];
    else
        rk2=[rk2;rkperN2'];
    end
    rkiter2=[rkiter2;iter];

    figure;
    plot(x,uexact,'-b',x,[0;A3\b3;0],'-.g',x,u1,'-.r');
    xlabel('x');
    ylabel('y');
    title(['Solution to the PDE at iteration k=' num2str(iter)...
        ' for N=' num2str(N)]);
    legend('exact','discretized','Biorthogonal CG','Location','SouthWest');
end

slope=polyfit(log10(Nvalues),log10(eN),1);
figure
plot(log10(Nvalues),log10(eN),'-.or')
xlabel('log10 ( N )')
ylabel('log10 ( max ( e^N ) )')
title(['Discretization error of the PDE for A3, slope = ' num2str(slope(1))])

figure
plot(1:rkiter1(1),log10(rk1(1,1:rkiter1(1))),1:rkiter1(2),...
    log10(rk1(2,1:rkiter1(2))),1:rkiter1(3),log10(rk1(3,1:rkiter1(3))),...
    1:rkiter1(4),log10(rk1(4,1:rkiter1(4))));
xlabel('k');
ylabel('log10 ( norm ( r_k ) )');
title('L2 norm of the Biorthogonal CG residual 1 value at each iteration');
legend(num2str(Nvalues(1)),num2str(Nvalues(2)),num2str(Nvalues(3)),...
    num2str(Nvalues(4)),'Location','SouthWest');
figure
plot(1:rkiter2(1),log10(rk2(1,1:rkiter2(1))),1:rkiter2(2),...
    log10(rk2(2,1:rkiter2(2))),1:rkiter2(3),log10(rk2(3,1:rkiter2(3))),...
    1:rkiter2(4),log10(rk2(4,1:rkiter2(4))));
xlabel('k');
ylabel('log10 ( norm ( s_k ) )');
title('L2 norm of the Biorthogonal CG residual 2 at each iteration');
legend(num2str(Nvalues(1)),num2str(Nvalues(2)),num2str(Nvalues(3)),...
    num2str(Nvalues(4)),'Location','SouthWest');
