% yasser@ucla.edu math 270c hw 3 Incomplete Cholesky (lower) Spring 2011
%
function L=mycholl(A)
L=tril(A);
n=size(A,1);
for k=1:n
    L(k,k)=sqrt(L(k,k));
    for i=k+1:n
        if A(i,k) ~= 0
            L(i,k)=L(i,k)/L(k,k);
        end
    end
    for j=k+1:n
        for i=j:n
            if A(i,j) ~= 0
                L(i,j)=L(i,j)-L(i,k)*L(j,k);
            end
        end
    end
end
