% yasser@ucla.edu math 270c hw 2 GMRES with restarts Spring 2011
%

close all
clear all
format longG
tol=1e-10;
maxiter=2;
N=3;
    
    x=[1 2 3];
    A3=[1 1 1;0 1 3;0 0 1];
    b3=[2;-4;1];
    
    [u3 flag relres iter r]=mygmres_restart(A3,b3,1,tol,maxiter);
    %[u3 flag relres iter r]=gmres(A3,b3,2,tol,maxiter);
     
    if maxiter > size(r,1)
        fill=zeros(1,maxiter-size(r,1));
        fill(1,:)=r(end);
        r=[r;fill];
    else
        r=r(1:maxiter);
    end
    figure;
    plot(x,A3\b3,'-.og',x,u3,'-.or');
    xlabel('x');
    ylabel('y');
    title(['Solution to Ax=b at iteration k=' num2str(iter(1,2))...
        ' for N=' num2str(N)]);
    legend('Exact','GMRES with restarts','Location','SouthWest');


figure;
plot(1:iter(1,1),log10(r(1:iter(1,1))),'-.o');
xlabel('k');
ylabel('log10 ( norm ( r_k ) )');
title('L2 norm of the GMRES with restarts residual at each iteration');
legend(num2str(N),'Location','SouthWest');
