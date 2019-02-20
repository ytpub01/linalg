% yasser@ucla.edu math 270c hw 3 Incomplete Cholesky (upper) Spring 2011
%
function R=mycholu(A)	%A spd matrix
R=triu(A);
n=size(A,1);
for k=1:n
    for j=k+1:n
        R(j,j:n)=R(j,j:n)-R(k,j:n)*R(k,j)/R(k,k);
    end
    R(k,k:n)=R(k,k:n)/sqrt(R(k,k));
end
