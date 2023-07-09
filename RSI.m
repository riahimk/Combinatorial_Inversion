function [L,U,iL,iU,iA]=RSI(A)
% Description: This function computes the UL factorization for a given square matrix A of size nxn. It also computes the LU factorization for the inverse of A.
% Inputs:
%   - A: The input square matrix of size nxn.
% Outputs:
%   - L: Lower triangular matrix of the UL factorization of A.
%   - U: Upper triangular matrix of the UL factorization of A.
%   - iL: Lower triangular matrix of the LU factorization of the inverse of A.
%   - iU: Upper triangular matrix of the LU factorization of the inverse of A.
% Author: Mohamed Kamel Riahi
% Affiliation: Mathematics Department, Khalifa University
% Date: [01/02/2023]



if nargin <1
    A=rand(6);
end

n  = size(A,1);

I  = eye(n);
Up = I;

[Ls,Us]=split(A);
for i=1:n
    Up=Up*Us;
    iU=CRIT(Us);
    [Ls,Us]=split(I+iU*Ls);
end
U  = Up;
L  = I+iU*Ls;
iU = CRIT(Up);
iL = CRIT(L')';
iA = iL*iU;

    function [Ls,Us]=split(B)
        Us = triu(B);
        Ls = tril(B,-1);
    end

end