function [S]=CRITe(R)
n=size(R,1);
S=eye(n,n);
for i=1:n-1
    S(i,i+1) = - R(i,i+1);
    for j=i+2:n
       S(i,j) = - S(i,i+1:j-1)*R(i+1:j-1,j) - R(i,j);
    end
end
end
