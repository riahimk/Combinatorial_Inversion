function [S,K,U,L] = SKUL(A)

if nargin < 1
    A=rand(5);A=A'*A +eye(size(A));
end

[n,~] = size(A);
L=zeros(n);
U=eye(n);
S=eye(n);
K=eye(n);



L(:, 1) = A(:, 1);

U(1, 2:n) = A(1, 2:n) / L(1, 1);

S(1,2)=-U(1,2);

K(1,1) = 1/L(1,1);

for i =2:n
    for j =2:i
        L(i, j) = A(i, j) - L(i, 1:j - 1) * U(1:j - 1, j);
    end
    
    K(i,1:i-1)= - L(i,1:i)*K(1:i,1:i-1)/L(i,i);
    
    K(:,i)  = K(:,i)/L(i,i);

    U(i, i+1:n) = (A(i, i+1:n) - L(i, 1:i - 1) * U(1:i - 1, i+1:n)) / L(i, i);
    
    if i<n
        for l=1:i
            S(l,i+1)= - S(l,l:i)*U(l:i,i+1);
        end
    end
end

end
